/**
 * ML-MAWS: Maximum Likelihood Phylogeny Estimation Using
 *          Minimal Absent Word Characters
 *
 * main.cpp - Main pipeline orchestrator
 *
 * Usage:
 *   ml-maws -i <input_dir_or_fasta> -o <output_dir> [options]
 *
 */

#include "FastaReader.h"
#include "MAWExtractor.h"
#include "MatrixBuilder.h"
#include "EntropySelector.h"
#include "ComplexityTracker.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <filesystem>
#include <cstdlib>
#include <sstream>
#include <algorithm>

namespace fs = std::filesystem;
using namespace mlmaws;

// ========== Configuration ==========
struct Config
{
    std::string inputPath;
    std::string outputDir = "./ml_maws_output";
    int fixedLength = 0;
    int lmin = 3;
    int lmax = 12;
    bool strandAware = false;
    bool skipML = false;
    bool useIQTree = false;
    bool adaptiveRange = true; // auto-compute lmin/lmax from sequence lengths
    bool multiLength = true;   // combine MAWs from top-K entropy lengths
    int topK = 3;              // number of lengths for multi-length mode
    bool fastML = false;       // fast ML mode for large datasets (>50 species)
    int threads = 4;
    int maxChars = 50000; // max characters in final matrix (0 = no limit)
};

void printUsage(const char *prog)
{
    std::cerr << "\nML-MAWS: Maximum Likelihood Phylogeny from Minimal Absent Words\n"
              << "=========================================================\n\n"
              << "Usage: " << prog << " -i <input> -o <output_dir> [options]\n\n"
              << "Required:\n"
              << "  -i <path>       Input directory of FASTA files or single multi-FASTA\n"
              << "  -o <path>       Output directory\n\n"
              << "Optional:\n"
              << "  -l <int>        Fixed MAW length (0 = auto-select, default: 0)\n"
              << "  --lmin <int>    Min MAW length for entropy search (default: 3)\n"
              << "  --lmax <int>    Max MAW length for entropy search (default: 12)\n"
              << "  --strand        Enable strand-aware MAW filtering\n"
              << "  --no-adaptive   Disable adaptive length range\n"
              << "  --no-multi      Use single best length only (no multi-length)\n"
              << "  --topk <int>    Lengths for multi-length aggregation (default: 3)\n"
              << "  --fast           Fast ML mode (for >50 species datasets)\n"
              << "  --no-raxml      Skip ML tree estimation (only produce PHYLIP matrix)\n"
              << "  --iqtree        Use IQ-TREE instead of RAxML\n"
              << "  --max-chars <n> Max characters in matrix (default: 50000, 0=unlimited)\n"
              << "  --threads <n>   Number of threads for ML tool (default: 4)\n"
              << "  -h, --help      Show this help\n\n";
}

Config parseArgs(int argc, char *argv[])
{
    Config cfg;
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if ((arg == "-i") && i + 1 < argc)
            cfg.inputPath = argv[++i];
        else if ((arg == "-o") && i + 1 < argc)
            cfg.outputDir = argv[++i];
        else if ((arg == "-l") && i + 1 < argc)
            cfg.fixedLength = std::atoi(argv[++i]);
        else if ((arg == "--lmin") && i + 1 < argc)
            cfg.lmin = std::atoi(argv[++i]);
        else if ((arg == "--lmax") && i + 1 < argc)
            cfg.lmax = std::atoi(argv[++i]);
        else if (arg == "--strand")
            cfg.strandAware = true;
        else if (arg == "--no-adaptive")
            cfg.adaptiveRange = false;
        else if (arg == "--no-multi")
            cfg.multiLength = false;
        else if ((arg == "--topk") && i + 1 < argc)
            cfg.topK = std::atoi(argv[++i]);
        else if (arg == "--fast")
            cfg.fastML = true;
        else if (arg == "--no-raxml")
            cfg.skipML = true;
        else if (arg == "--iqtree")
            cfg.useIQTree = true;
        else if ((arg == "--max-chars") && i + 1 < argc)
            cfg.maxChars = std::atoi(argv[++i]);
        else if ((arg == "--threads") && i + 1 < argc)
            cfg.threads = std::atoi(argv[++i]);
        else if (arg == "-h" || arg == "--help")
        {
            printUsage(argv[0]);
            exit(0);
        }
        else
        {
            std::cerr << "[WARN] Unknown argument: " << arg << std::endl;
        }
    }
    return cfg;
}

// ========== Helper: check if file exists ==========
static bool fileExists(const std::string &path)
{
    return fs::exists(path);
}

// ========== Helper: remove file if exists ==========
static void removeFileIfExists(const std::string &path)
{
    if (fs::exists(path))
        fs::remove(path);
}

// ========== Helper: copy file ==========
static void copyFile(const std::string &src, const std::string &dst)
{
    fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
}

// ========== Helper: get absolute path ==========
static std::string getAbsolutePath(const std::string &path)
{
    return fs::absolute(path).string();
}

// ========== Main Pipeline ==========
int main(int argc, char *argv[])
{
    auto totalStart = std::chrono::high_resolution_clock::now();

    // Complexity tracker for time & space analysis
    ComplexityTracker tracker;

    std::cerr << "\n"
              << "============================================================\n"
              << "  ML-MAWS: Maximum Likelihood Phylogeny Estimation\n"
              << "  Using Minimal Absent Word Characters\n"
              << "============================================================\n\n";

    if (argc < 3)
    {
        printUsage(argv[0]);
        return 1;
    }
    Config cfg = parseArgs(argc, argv);

    if (cfg.inputPath.empty())
    {
        std::cerr << "[ERROR] No input path specified. Use -i <path>." << std::endl;
        return 1;
    }

    fs::create_directories(cfg.outputDir);

    // ============================================================
    // STEP 1: Read input sequences
    // ============================================================
    std::cerr << "[STEP 1/5] Reading input sequences..." << std::endl;
    tracker.startStep("Read FASTA", "O(m * n)", "");

    std::vector<FastaRecord> records;
    if (fs::is_directory(cfg.inputPath))
    {
        records = FastaReader::readDirectory(cfg.inputPath);
    }
    else
    {
        records = FastaReader::readFile(cfg.inputPath);
    }

    if (records.size() < 3)
    {
        std::cerr << "[ERROR] Need at least 3 species. Found: " << records.size() << std::endl;
        return 1;
    }

    // Compute average sequence length for complexity reporting
    long long totalBP = 0;
    for (const auto &r : records)
        totalBP += r.sequence.size();
    int avgSeqLen = (int)(totalBP / records.size());

    tracker.endStep();
    std::cerr << "[STEP 1] Done. " << records.size() << " species (" << tracker.getTotalTimeSeconds() << "s)\n";

    for (size_t i = 0; i < records.size(); i++)
    {
        std::cerr << "  [" << (i + 1) << "] " << records[i].name
                  << " (" << records[i].sequence.size() << " bp)\n";
    }
    std::cerr << std::endl;

    // ============================================================
    // FIX 1: Adaptive length range (solves short sequences)
    // ============================================================
    int m = (int)records.size();

    if (cfg.adaptiveRange && cfg.fixedLength == 0)
    {
        EntropySelector::computeAdaptiveRange(avgSeqLen, m, cfg.lmin, cfg.lmax);
    }

    // Auto-enable fast ML for large datasets (>50 species)
    if (m > 50 && !cfg.fastML)
    {
        std::cerr << "[INFO] Auto-enabling fast ML mode for " << m << " species\n";
        cfg.fastML = true;
    }

    // ============================================================
    // STEP 2: Extract MAWs (using adaptive range)
    // ============================================================
    std::cerr << "[STEP 2/5] Extracting Minimal Absent Words (range=["
              << cfg.lmin << "," << cfg.lmax << "])..." << std::endl;
    std::string step2params = "m=" + std::to_string(m) + ", n_avg=" + std::to_string(avgSeqLen) + ", sigma=4";
    tracker.startStep("MAW Extraction", "O(m * n * sigma)", step2params);

    MAWExtractor extractor(cfg.lmin, cfg.lmax);
    std::vector<std::string> speciesNames(m);
    std::vector<std::vector<std::string>> allMAWSets(m);

    for (int i = 0; i < m; i++)
    {
        speciesNames[i] = records[i].name;
        auto t0 = std::chrono::high_resolution_clock::now();
        if (cfg.strandAware)
            allMAWSets[i] = extractor.extractStrandAware(records[i].sequence);
        else
            allMAWSets[i] = extractor.extract(records[i].sequence);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::cerr << "  " << speciesNames[i] << ": " << allMAWSets[i].size()
                  << " MAWs (" << std::chrono::duration<double>(t1 - t0).count() << "s)\n";
    }

    tracker.endStep();
    std::cerr << "[STEP 2] Done (" << tracker.getTotalTimeSeconds() << "s)\n\n";

    // ============================================================
    // STEP 3: Optimal length selection (entropy-based)
    // FIX 2: Multi-length aggregation (solves similar species)
    // ============================================================
    std::vector<int> selectedLengths;

    if (cfg.fixedLength != 0)
    {
        selectedLengths.push_back(cfg.fixedLength);
        std::cerr << "[STEP 3/5] Using fixed MAW length: " << cfg.fixedLength << "\n\n";
    }
    else
    {
        std::cerr << "[STEP 3/5] Selecting optimal MAW length(s) via entropy..." << std::endl;
        std::string step3params = "L_range=[" + std::to_string(cfg.lmin) + ".." + std::to_string(cfg.lmax) + "]";
        tracker.startStep("Entropy Selection", "O(L * m * |U_max|)", step3params);

        std::vector<EntropyResult> entropyResults;

        for (int len = cfg.lmin; len <= cfg.lmax; len++)
        {
            std::vector<std::vector<std::string>> lengthMAWSets(m);
            for (int i = 0; i < m; i++)
            {
                for (const auto &w : allMAWSets[i])
                {
                    if ((int)w.size() == len)
                        lengthMAWSets[i].push_back(w);
                }
            }

            MatrixBuilder tempBuilder;
            tempBuilder.build(speciesNames, lengthMAWSets, true);

            if (tempBuilder.numCharacters() == 0)
            {
                entropyResults.push_back({len, 0.0, 0});
                continue;
            }

            double entropy = EntropySelector::computeEntropy(tempBuilder.getMatrix(), m);
            entropyResults.push_back({len, entropy, tempBuilder.numCharacters()});
        }

        // Multi-length or single-length selection
        if (cfg.multiLength)
        {
            // FIX 2: Combine MAWs from multiple high-entropy lengths
            selectedLengths = EntropySelector::selectTopLengths(entropyResults, cfg.topK, 5);
            if (selectedLengths.empty())
            {
                // Fallback to single best
                int best = EntropySelector::selectBestLength(entropyResults);
                if (best > 0)
                    selectedLengths.push_back(best);
            }
        }
        else
        {
            int best = EntropySelector::selectBestLength(entropyResults);
            if (best > 0)
                selectedLengths.push_back(best);
        }

        // Save entropy results
        std::string entropyFile = cfg.outputDir + "/entropy_results.tsv";
        std::ofstream efile(entropyFile);
        efile << "Length\tCharacters\tEntropy\n";
        for (const auto &r : entropyResults)
            efile << r.length << "\t" << r.numCharacters << "\t" << r.entropy << "\n";
        efile.close();

        tracker.endStep();
        std::cerr << "[STEP 3] Done\n\n";
    }

    if (selectedLengths.empty())
    {
        std::cerr << "[ERROR] Could not determine optimal MAW length." << std::endl;
        return 1;
    }

    // ============================================================
    // STEP 4: Build final binary character matrix
    // FIX 2 continued: Aggregate MAWs from all selected lengths
    // ============================================================
    std::string lengthsStr;
    for (size_t i = 0; i < selectedLengths.size(); i++)
    {
        if (i > 0)
            lengthsStr += ",";
        lengthsStr += std::to_string(selectedLengths[i]);
    }
    std::cerr << "[STEP 4/5] Building binary matrix (lengths={" << lengthsStr << "})..." << std::endl;
    tracker.startStep("Matrix Construction", "O(m * |U|)", "lengths={" + lengthsStr + "}");

    // Collect MAWs from ALL selected lengths
    std::vector<std::vector<std::string>> filteredMAWSets(m);
    for (int i = 0; i < m; i++)
    {
        for (const auto &w : allMAWSets[i])
        {
            int wlen = (int)w.size();
            for (int sl : selectedLengths)
            {
                if (wlen == sl)
                {
                    filteredMAWSets[i].push_back(w);
                    break;
                }
            }
        }
        // Re-sort after aggregating multiple lengths
        std::sort(filteredMAWSets[i].begin(), filteredMAWSets[i].end());
    }
    allMAWSets.clear();
    allMAWSets.shrink_to_fit();

    MatrixBuilder builder;
    builder.build(speciesNames, filteredMAWSets, true);

    if (builder.numCharacters() == 0)
    {
        std::cerr << "[ERROR] No variable characters. Try different length range." << std::endl;
        return 1;
    }

    // Character capping: keep most informative columns if too many
    if (cfg.maxChars > 0 && builder.numCharacters() > cfg.maxChars)
    {
        std::cerr << "[INFO] Matrix has " << builder.numCharacters()
                  << " characters, capping to " << cfg.maxChars
                  << " most informative..." << std::endl;
        builder.capCharacters(cfg.maxChars);
    }

    std::cerr << "[INFO] Matrix: " << m << " species x " << builder.numCharacters()
              << " characters (" << selectedLengths.size() << " length(s))" << std::endl;

    std::string phylipFile = cfg.outputDir + "/ml_maws_matrix.phy";
    builder.writePhylip(phylipFile);

    std::string tsvFile = cfg.outputDir + "/ml_maws_matrix.tsv";
    builder.writeTSV(tsvFile);

    tracker.endStep();
    std::cerr << "[STEP 4] Done\n\n";

    // ============================================================
    // STEP 5: Run ML tree estimation
    // FIX 3: Fast ML mode for large datasets
    // ============================================================
    if (!cfg.skipML)
    {
        std::cerr << "[STEP 5/5] Running ML tree estimation";
        if (cfg.fastML)
            std::cerr << " (FAST mode)";
        std::cerr << "..." << std::endl;
        std::string step5params = "m=" + std::to_string(m) + ", |U|=" + std::to_string(builder.numCharacters());
        if (cfg.fastML)
            step5params += ", mode=fast";
        tracker.startStep("ML Inference", "NP-hard (heuristic)", step5params);

        std::string treeFile;
        std::string cmdStr;
        int ret;

        if (cfg.useIQTree)
        {
            std::string prefix = cfg.outputDir + "/ml_maws_iqtree";
            if (cfg.fastML)
            {
                // Fast mode: fewer bootstrap, fast search
                cmdStr = "iqtree2 -s " + phylipFile +
                         " -st BIN -m MFP+ASC" +
                         " -bb 1000 -nt " + std::to_string(cfg.threads) +
                         " --prefix " + prefix + " -redo -fast";
            }
            else
            {
                cmdStr = "iqtree2 -s " + phylipFile +
                         " -st BIN -m MFP+ASC" +
                         " -bb 1000 -nt " + std::to_string(cfg.threads) +
                         " --prefix " + prefix + " -redo";
            }
            treeFile = prefix + ".treefile";
        }
        else
        {
            std::string runName = "ML_MAWS";
            std::vector<std::string> suffixes = {
                "bestTree", "bipartitions", "bipartitionsBranchLabels",
                "bootstrap", "info", "log", "parsimonyTree", "result"};
            for (const auto &suf : suffixes)
                removeFileIfExists(cfg.outputDir + "/RAxML_" + suf + "." + runName);

            std::string absOutDir = getAbsolutePath(cfg.outputDir);
            int bootstrapReps = cfg.fastML ? 20 : 100; // Fewer reps in fast mode

            cmdStr = "raxmlHPC-PTHREADS -f a -m BINGAMMA"
                     " -p 12345 -x 12345 -N " +
                     std::to_string(bootstrapReps) +
                     " -T " + std::to_string(cfg.threads) +
                     " -s " + phylipFile +
                     " -w " + absOutDir +
                     " -n " + runName;
            treeFile = cfg.outputDir + "/RAxML_bestTree." + runName;
        }

        std::cerr << "[CMD] " << cmdStr << std::endl;
        ret = std::system(cmdStr.c_str());

        if (ret != 0 && !cfg.useIQTree)
        {
            std::cerr << "[WARN] Threaded RAxML failed, trying serial..." << std::endl;
            std::string runName = "ML_MAWS";
            std::string absOutDir = getAbsolutePath(cfg.outputDir);
            cmdStr = "raxmlHPC -f a -m BINGAMMA -p 12345 -x 12345 -N 100"
                     " -s " +
                     phylipFile + " -w " + absOutDir + " -n " + runName;
            std::cerr << "[CMD] " << cmdStr << std::endl;
            ret = std::system(cmdStr.c_str());
        }

        tracker.endStep();

        if (ret == 0 && fileExists(treeFile))
        {
            std::string finalTree = cfg.outputDir + "/ML_MAWS_tree.newick";
            copyFile(treeFile, finalTree);
            std::cerr << "[STEP 5] Done\n"
                      << "[OUTPUT] Tree: " << finalTree << std::endl;
        }
        else
        {
            std::cerr << "[STEP 5] ML tool failed.\n"
                      << "[INFO] Matrix at: " << phylipFile << "\n"
                      << "[INFO] Run manually:\n"
                      << "  raxmlHPC -m BINGAMMA -p 12345 -s " << phylipFile << " -n T1\n"
                      << "  iqtree2 -s " << phylipFile << " -st BIN -m MFP+ASC -bb 1000\n";
        }
    }
    else
    {
        std::cerr << "[STEP 5/5] Skipped (--no-raxml).\n"
                  << "[INFO] Matrix: " << phylipFile << "\n";
    }

    // ============================================================
    // COMPLEXITY ANALYSIS OUTPUT
    // ============================================================

    // Print per-step time & memory summary table
    tracker.printSummary(m, avgSeqLen, builder.numCharacters());

    // Write detailed complexity report JSON
    std::string complexityFile = cfg.outputDir + "/complexity_report.json";
    tracker.writeJSON(complexityFile, m, avgSeqLen,
                      builder.numCharacters(), selectedLengths[0], cfg.strandAware);
    std::cerr << "  Complexity report: " << complexityFile << "\n";

    // Final summary
    auto totalEnd = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration<double>(totalEnd - totalStart).count();

    std::cerr << "\n============================================================\n"
              << "  SUMMARY\n"
              << "============================================================\n"
              << "  Species:      " << m << "\n"
              << "  MAW length:   {" << lengthsStr << "}\n"
              << "  Characters:   " << builder.numCharacters() << "\n"
              << "  Strand-aware: " << (cfg.strandAware ? "Yes" : "No") << "\n"
              << "  Total time:   " << totalTime << " sec\n"
              << "  Peak memory:  " << (tracker.getPeakMemoryKB() / 1024.0) << " MB\n"
              << "  Output dir:   " << cfg.outputDir << "\n"
              << "============================================================\n\n";

    return 0;
}
