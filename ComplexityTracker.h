/**
 * ComplexityTracker.h - Time and Space Complexity Measurement
 *
 * Tracks per-step wall-clock time, peak memory (RSS), and
 * reports theoretical vs. empirical complexity.
 *
 * Works on Linux (reads /proc/self/status) and falls back
 * to getrusage() on other POSIX systems.
 */

#ifndef COMPLEXITY_TRACKER_H
#define COMPLEXITY_TRACKER_H

#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef _WIN32
    #include <windows.h>
    #include <psapi.h>
#else
    #include <sys/resource.h>
    #include <unistd.h>
#endif

namespace mlmaws {

/**
 * A single step's performance record.
 */
struct StepRecord {
    std::string name;
    double time_seconds;
    long peak_memory_kb;      // peak RSS in KB at end of step
    long memory_delta_kb;     // memory increase during this step
    std::string complexity;   // theoretical complexity string
    std::string parameters;   // e.g., "m=25, n=16000, sigma=4"
};

/**
 * Tracks time and memory usage across pipeline steps.
 * Outputs a JSON summary at the end.
 */
class ComplexityTracker {
public:
    ComplexityTracker() : initialMemory_(getCurrentMemoryKB()) {}

    /** Start timing a new step. */
    void startStep(const std::string& name, const std::string& complexity = "",
                   const std::string& parameters = "") {
        currentStep_.name = name;
        currentStep_.complexity = complexity;
        currentStep_.parameters = parameters;
        currentStep_.peak_memory_kb = 0;
        currentStep_.memory_delta_kb = 0;
        stepStartMemory_ = getCurrentMemoryKB();
        stepStartTime_ = std::chrono::high_resolution_clock::now();
    }

    /** End the current step and record results. */
    void endStep() {
        auto endTime = std::chrono::high_resolution_clock::now();
        currentStep_.time_seconds = std::chrono::duration<double>(
            endTime - stepStartTime_).count();
        long endMemory = getCurrentMemoryKB();
        currentStep_.peak_memory_kb = endMemory;
        currentStep_.memory_delta_kb = endMemory - stepStartMemory_;
        steps_.push_back(currentStep_);
    }

    /** Get current memory usage (peak RSS) in KB. */
    static long getCurrentMemoryKB() {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
            return (long)(pmc.PeakWorkingSetSize / 1024);
        }
        return 0;
#else
        // Try /proc/self/status first (Linux)
        std::ifstream status("/proc/self/status");
        if (status.is_open()) {
            std::string line;
            while (std::getline(status, line)) {
                if (line.substr(0, 6) == "VmRSS:") {
                    std::istringstream iss(line.substr(6));
                    long val;
                    iss >> val;
                    return val; // already in KB
                }
            }
        }
        // Fallback: getrusage
        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) == 0) {
            #ifdef __APPLE__
                return usage.ru_maxrss / 1024; // macOS reports bytes
            #else
                return usage.ru_maxrss; // Linux reports KB
            #endif
        }
        return 0;
#endif
    }

    /** Get peak memory across all steps. */
    long getPeakMemoryKB() const {
        long peak = 0;
        for (const auto& s : steps_) {
            peak = std::max(peak, s.peak_memory_kb);
        }
        return peak;
    }

    /** Get total time across all steps. */
    double getTotalTimeSeconds() const {
        double total = 0;
        for (const auto& s : steps_) {
            total += s.time_seconds;
        }
        return total;
    }

    /** Print summary table to stderr. */
    void printSummary(int numSpecies, int seqLength, int numCharacters) const {
        std::cerr << "\n============================================================\n"
                  << "  TIME AND SPACE COMPLEXITY ANALYSIS\n"
                  << "============================================================\n\n";

        // Parameters
        std::cerr << "  Input Parameters:\n"
                  << "    m (species):     " << numSpecies << "\n"
                  << "    n (avg seq len): " << seqLength << "\n"
                  << "    |U| (characters):" << numCharacters << "\n"
                  << "    sigma (alphabet):" << 4 << "\n\n";

        // Per-step table
        std::cerr << "  " << std::left
                  << std::setw(30) << "Step"
                  << std::setw(12) << "Time (s)"
                  << std::setw(14) << "Memory (MB)"
                  << std::setw(12) << "Delta (MB)"
                  << "Complexity\n";
        std::cerr << "  " << std::string(90, '-') << "\n";

        for (const auto& s : steps_) {
            std::cerr << "  " << std::left
                      << std::setw(30) << s.name
                      << std::setw(12) << std::fixed << std::setprecision(3) << s.time_seconds
                      << std::setw(14) << std::fixed << std::setprecision(2)
                      << (s.peak_memory_kb / 1024.0)
                      << std::setw(12) << std::fixed << std::setprecision(2)
                      << (s.memory_delta_kb / 1024.0)
                      << s.complexity << "\n";
            if (!s.parameters.empty()) {
                std::cerr << "  " << std::setw(30) << "" << "  (" << s.parameters << ")\n";
            }
        }

        std::cerr << "  " << std::string(90, '-') << "\n";
        std::cerr << "  " << std::left
                  << std::setw(30) << "TOTAL"
                  << std::setw(12) << std::fixed << std::setprecision(3) << getTotalTimeSeconds()
                  << std::setw(14) << std::fixed << std::setprecision(2)
                  << (getPeakMemoryKB() / 1024.0)
                  << "\n\n";

        // Theoretical complexity summary
        std::cerr << "  Theoretical Complexity:\n"
                  << "    Step 1 (Read FASTA):     O(m * n)\n"
                  << "    Step 2 (Extract MAWs):   O(m * n * sigma)\n"
                  << "    Step 3 (Entropy select):  O(L_range * m * |U_max|)\n"
                  << "    Step 4 (Build matrix):   O(m * |U|)\n"
                  << "    Step 5 (ML inference):   NP-hard (RAxML/IQ-TREE heuristic)\n"
                  << "    Overall:                 O(m * n * sigma) + ML(m, |U|)\n\n";
    }

    /** Write full results to JSON file. */
    void writeJSON(const std::string& filepath, int numSpecies, int seqLength,
                   int numCharacters, int mawLength, bool strandAware) const {
        std::ofstream out(filepath);
        if (!out.is_open()) return;

        out << "{\n";
        out << "  \"parameters\": {\n"
            << "    \"num_species\": " << numSpecies << ",\n"
            << "    \"avg_sequence_length\": " << seqLength << ",\n"
            << "    \"num_characters\": " << numCharacters << ",\n"
            << "    \"maw_length\": " << mawLength << ",\n"
            << "    \"strand_aware\": " << (strandAware ? "true" : "false") << ",\n"
            << "    \"alphabet_size\": 4\n"
            << "  },\n";

        out << "  \"steps\": [\n";
        for (size_t i = 0; i < steps_.size(); i++) {
            const auto& s = steps_[i];
            out << "    {\n"
                << "      \"name\": \"" << s.name << "\",\n"
                << "      \"time_seconds\": " << std::fixed << std::setprecision(6) << s.time_seconds << ",\n"
                << "      \"peak_memory_kb\": " << s.peak_memory_kb << ",\n"
                << "      \"memory_delta_kb\": " << s.memory_delta_kb << ",\n"
                << "      \"peak_memory_mb\": " << std::fixed << std::setprecision(2)
                << (s.peak_memory_kb / 1024.0) << ",\n"
                << "      \"complexity\": \"" << s.complexity << "\",\n"
                << "      \"parameters\": \"" << s.parameters << "\"\n"
                << "    }" << (i + 1 < steps_.size() ? "," : "") << "\n";
        }
        out << "  ],\n";

        out << "  \"totals\": {\n"
            << "    \"total_time_seconds\": " << std::fixed << std::setprecision(6) << getTotalTimeSeconds() << ",\n"
            << "    \"peak_memory_kb\": " << getPeakMemoryKB() << ",\n"
            << "    \"peak_memory_mb\": " << std::fixed << std::setprecision(2)
            << (getPeakMemoryKB() / 1024.0) << "\n"
            << "  },\n";

        out << "  \"theoretical_complexity\": {\n"
            << "    \"time\": \"O(m * n * sigma) + ML(m, |U|)\",\n"
            << "    \"space\": \"O(n * sigma) for SA + O(m * |U|) for matrix\",\n"
            << "    \"step1_read\": \"O(m * n)\",\n"
            << "    \"step2_maw_extract\": \"O(m * n * sigma)\",\n"
            << "    \"step3_entropy\": \"O(L_range * m * |U_max|)\",\n"
            << "    \"step4_matrix\": \"O(m * |U|)\",\n"
            << "    \"step5_ml\": \"NP-hard (heuristic)\"\n"
            << "  }\n";

        out << "}\n";
        out.close();
    }

private:
    std::vector<StepRecord> steps_;
    StepRecord currentStep_;
    long initialMemory_;
    long stepStartMemory_;
    std::chrono::high_resolution_clock::time_point stepStartTime_;
};

} // namespace mlmaws

#endif // COMPLEXITY_TRACKER_H
