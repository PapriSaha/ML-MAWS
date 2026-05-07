#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}/.."
TOOLS_DIR="${PROJECT_DIR}/tools"
THREADS=4

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_DIR="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        --ref) REF_TREE="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        *) shift ;;
    esac
done

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 -i <input_dir> -o <output_dir> --ref <ref.newick> [-t threads]"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Collection of FASTA files safely
shopt -s nullglob
FASTA_FILES=("$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fna)
shopt -u nullglob
NUM_SPECIES=${#FASTA_FILES[@]}

if [ "$NUM_SPECIES" -eq 0 ]; then
    echo "ERROR: No FASTA files (.fasta/.fa/.fna) found in $INPUT_DIR"
    echo "  Contents of $INPUT_DIR:"
    ls "$INPUT_DIR/" | head -10
    exit 1
fi

echo ""
echo "============================================================"
echo "  Running ALL Methods (${NUM_SPECIES} species)"
echo "  Input: $INPUT_DIR"
echo "  Output: $OUTPUT_DIR"
echo "  Reference: $REF_TREE"
echo "============================================================"
echo ""

# Results file
RESULTS_TSV="$OUTPUT_DIR/results.tsv"
echo -e "Method\tnRF\tnQD\tTime(s)\tMemory(MB)\tTree_File" > "$RESULTS_TSV"

# ============================================================
# Helper: compute nRF between two trees
# ============================================================
compute_nrf() {
    local ref_tree="$1"
    local est_tree="$2"

    if [ -z "$ref_tree" ] || [ ! -f "$ref_tree" ] || [ ! -f "$est_tree" ]; then
        echo "N/A"
        return
    fi

    python3 -c "
import dendropy
from dendropy.calculate import treecompare
tns = dendropy.TaxonNamespace()
ref = dendropy.Tree.get(path='$ref_tree', schema='newick', taxon_namespace=tns, rooting='default-unrooted')
est = dendropy.Tree.get(path='$est_tree', schema='newick', taxon_namespace=tns, rooting='default-unrooted')
rf = treecompare.symmetric_difference(ref, est)
n = len(tns)
max_rf = 2 * (n - 3)
print(f'{rf/max_rf:.6f}' if max_rf > 0 else '0.0')
" 2>/dev/null || echo "ERROR"
}

# ============================================================
# Helper: record result for a method
# ============================================================
record_result() {
    local method_name="$1"
    local tree_file="$2"
    local elapsed="$3"
    local peak_mem_mb="$4"

    local nrf
    nrf=$(compute_nrf "$REF_TREE" "$tree_file")

    echo "    Tree: $tree_file"
    echo "    Time: ${elapsed}s, Memory: ${peak_mem_mb} MB, nRF: $nrf"
    echo -e "${method_name}\t${nrf}\tN/A\t${elapsed}\t${peak_mem_mb}\t${tree_file}" >> "$RESULTS_TSV"
}

# ============================================================
# 1. ML-MAWS (use EXISTING trees from benchmark results)
# ============================================================
DS_DIR=$(dirname "$OUTPUT_DIR")  # e.g. benchmark/results/ecoli-shigella

# Check for existing IQ-TREE result
EXISTING_IQTREE="$DS_DIR/ml-maws-iqtree/ML_MAWS_tree.newick"
if [ -f "$EXISTING_IQTREE" ]; then
    echo "  [ML-MAWS-IQT] Using existing tree..."
    record_result "ML-MAWS-IQT" "$EXISTING_IQTREE" "pre-computed" "N/A"
else
    echo "  [ML-MAWS-IQT] No existing tree found at $EXISTING_IQTREE"
fi

# Check for existing nostrand result
EXISTING_NOSTRAND="$DS_DIR/ml-maws-nostrand/ML_MAWS_tree.newick"
if [ -f "$EXISTING_NOSTRAND" ]; then
    echo "  [ML-MAWS-NoStrand] Using existing tree..."
    record_result "ML-MAWS-NoStrand" "$EXISTING_NOSTRAND" "pre-computed" "N/A"
else
    echo "  [ML-MAWS-NoStrand] No existing tree found."
fi

# Check for existing RAxML result
EXISTING_RAXML="$DS_DIR/ml-maws-raxml/ML_MAWS_tree.newick"
if [ -f "$EXISTING_RAXML" ]; then
    echo "  [ML-MAWS-RAxML] Using existing tree..."
    record_result "ML-MAWS-RAxML" "$EXISTING_RAXML" "pre-computed" "N/A"
fi

# ============================================================
# 2. Mash + NJ
# ============================================================
if command -v mash &>/dev/null; then
    MASH_OUT="$OUTPUT_DIR/Mash"
    mkdir -p "$MASH_OUT"
    echo "  [Mash] Running..."

    START=$(date +%s)
    mash triangle -p "$THREADS" "${FASTA_FILES[@]}" > "$MASH_OUT/dist_raw.phylip" 2>"$MASH_OUT/stderr.log"
    # Strip full paths to basenames (remove path and .fasta extension)
    sed -E 's|/[^ ]*/([^/]+)\.fasta|\1|g; s|/[^ ]*/([^/]+)\.fa|\1|g; s|/[^ ]*/([^/]+)\.fna|\1|g' \
        "$MASH_OUT/dist_raw.phylip" > "$MASH_OUT/dist.phylip"
    python3 "$SCRIPT_DIR/dist_to_nj.py" "$MASH_OUT/dist.phylip" "$MASH_OUT/tree.newick" phylip
    END=$(date +%s)
    ELAPSED=$((END - START))

    record_result "Mash" "$MASH_OUT/tree.newick" "$ELAPSED" "N/A"
else
    echo "  [Mash] SKIPPED (not installed)"
fi

# ============================================================
# 3. Skmer + NJ
# ============================================================
if command -v skmer &>/dev/null; then
    SKMER_OUT="$OUTPUT_DIR/Skmer"
    mkdir -p "$SKMER_OUT"
    echo "  [Skmer] Running..."

    START=$(date +%s)
    skmer reference "$INPUT_DIR" -o "$SKMER_OUT" 2>"$SKMER_OUT/stderr.log" || true
    # Find the library directory skmer created
    SKMER_LIB=$(find "$SKMER_OUT" -name "CONFIG" -exec dirname {} \; 2>/dev/null | head -1)
    if [ -n "$SKMER_LIB" ]; then
        skmer distance "$SKMER_LIB" -o "$SKMER_OUT/dist.txt" 2>>"$SKMER_OUT/stderr.log" || true
    fi
    if [ -f "$SKMER_OUT/dist.txt" ]; then
        python3 "$SCRIPT_DIR/dist_to_nj.py" "$SKMER_OUT/dist.txt" "$SKMER_OUT/tree.newick" tsv
    fi
    END=$(date +%s)
    ELAPSED=$((END - START))

    record_result "Skmer" "$SKMER_OUT/tree.newick" "$ELAPSED" "N/A"
else
    echo "  [Skmer] SKIPPED (not installed)"
fi

# ============================================================
# 4. andi + NJ
# ============================================================
if command -v andi &>/dev/null; then
    ANDI_OUT="$OUTPUT_DIR/andi"
    mkdir -p "$ANDI_OUT"
    echo "  [andi] Running..."

    START=$(date +%s)
    # andi takes individual FASTA files as arguments (not concatenated)
    andi -j "${FASTA_FILES[@]}" > "$ANDI_OUT/dist.phylip" 2>"$ANDI_OUT/stderr.log"
    # Strip full paths to basenames in the phylip output
    sed -E 's|/[^ ]*/([^/]+)\.fasta|\1|g; s|/[^ ]*/([^/]+)\.fa|\1|g; s|/[^ ]*/([^/]+)\.fna|\1|g' \
        -i "$ANDI_OUT/dist.phylip"
    python3 "$SCRIPT_DIR/dist_to_nj.py" "$ANDI_OUT/dist.phylip" "$ANDI_OUT/tree.newick" phylip
    END=$(date +%s)
    ELAPSED=$((END - START))

    record_result "andi" "$ANDI_OUT/tree.newick" "$ELAPSED" "N/A"
else
    echo "  [andi] SKIPPED (not installed)"
fi

# ============================================================
# 5. CD-MAWS-SA + NJ
# ============================================================
CDMAWS_MAWGEN="${TOOLS_DIR}/cd-maws-sa/cd-maws-sa/maw-gen"
CDMAWS_DIST="${TOOLS_DIR}/cd-maws-sa/cd-maws-sa/dist-calc"
if [ -f "$CDMAWS_MAWGEN" ] && [ -f "$CDMAWS_DIST" ]; then
    CD_OUT="$OUTPUT_DIR/CD-MAWS"
    mkdir -p "$CD_OUT"
    echo "  [CD-MAWS] Running..."

    MULTI_FASTA="$CD_OUT/all.fasta"
    cat "${FASTA_FILES[@]}" > "$MULTI_FASTA"

    START=$(date +%s)
    "$CDMAWS_MAWGEN" "$MULTI_FASTA" "$CD_OUT/maws.txt" 2 12 2>"$CD_OUT/stderr.log" || true
    "$CDMAWS_DIST" "$CD_OUT/maws.txt" "$CD_OUT/dist.csv" 2>>"$CD_OUT/stderr.log" || true
    if [ -f "$CD_OUT/dist.csv" ]; then
        python3 "$SCRIPT_DIR/dist_to_nj.py" "$CD_OUT/dist.csv" "$CD_OUT/tree.newick" csv
    fi
    END=$(date +%s)
    ELAPSED=$((END - START))

    record_result "CD-MAWS" "$CD_OUT/tree.newick" "$ELAPSED" "N/A"
else
    echo "  [CD-MAWS] SKIPPED (not compiled at $CDMAWS_MAWGEN)"
fi

# ============================================================
# 6. FFP + NJ
# ============================================================
if command -v ffpry &>/dev/null; then
    FFP_OUT="$OUTPUT_DIR/FFP"
    mkdir -p "$FFP_OUT"
    echo "  [FFP] Running..."

    START=$(date +%s)
    ffpry -l 18 "${FASTA_FILES[@]}" | ffpcol | ffprwn | ffpjsd > "$FFP_OUT/dist.phylip" 2>"$FFP_OUT/stderr.log"
    python3 "$SCRIPT_DIR/dist_to_nj.py" "$FFP_OUT/dist.phylip" "$FFP_OUT/tree.newick" phylip
    END=$(date +%s)
    ELAPSED=$((END - START))

    record_result "FFP" "$FFP_OUT/tree.newick" "$ELAPSED" "N/A"
else
    echo "  [FFP] SKIPPED (not installed)"
fi

# ============================================================
# 7. Co-phylog + NJ
# ============================================================
COPHYLOG="${TOOLS_DIR}/co-phylog/co-phylog"
if [ -f "$COPHYLOG" ]; then
    COPHY_OUT="$OUTPUT_DIR/Co-phylog"
    mkdir -p "$COPHY_OUT"
    echo "  [Co-phylog] Running..."

    MULTI_FASTA="$COPHY_OUT/all.fasta"
    cat "${FASTA_FILES[@]}" > "$MULTI_FASTA"

    START=$(date +%s)
    "$COPHYLOG" "$MULTI_FASTA" > "$COPHY_OUT/dist.phylip" 2>"$COPHY_OUT/stderr.log" || true
    if [ -f "$COPHY_OUT/dist.phylip" ]; then
        python3 "$SCRIPT_DIR/dist_to_nj.py" "$COPHY_OUT/dist.phylip" "$COPHY_OUT/tree.newick" phylip
    fi
    END=$(date +%s)
    ELAPSED=$((END - START))

    record_result "Co-phylog" "$COPHY_OUT/tree.newick" "$ELAPSED" "N/A"
else
    echo "  [Co-phylog] SKIPPED (not compiled)"
fi

# ============================================================
# 8. kSNP4 (or kSNP3)
# ============================================================
KSNP_BIN=""
if command -v kSNP4 &>/dev/null; then
    KSNP_BIN="kSNP4"
elif command -v kSNP3 &>/dev/null; then
    KSNP_BIN="kSNP3"
fi

if [ -n "$KSNP_BIN" ]; then
    KSNP_OUT="$OUTPUT_DIR/kSNP"
    mkdir -p "$KSNP_OUT"
    echo "  [kSNP] Running with $KSNP_BIN..."

    # Create genomes.list (tab-separated: path\tname)
    GENOME_LIST="$KSNP_OUT/genomes.list"
    > "$GENOME_LIST"
    for f in "${FASTA_FILES[@]}"; do
        base=$(basename "$f" | sed 's/\.\(fasta\|fa\|fna\)$//')
        echo -e "$(realpath "$f")\t${base}" >> "$GENOME_LIST"
    done

    START=$(date +%s)
    # Determine optimal k
    KCHOOSER="Kchooser4"
    command -v "$KCHOOSER" &>/dev/null || KCHOOSER="Kchooser"
    if command -v "$KCHOOSER" &>/dev/null; then
        "$KCHOOSER" -in "$GENOME_LIST" 2>"$KSNP_OUT/kchooser_stderr.log" || true
        OPT_K=$(grep -oP 'The optimum value of K is \K[0-9]+' Kchooser.report 2>/dev/null || echo "21")
        mv Kchooser.report "$KSNP_OUT/" 2>/dev/null || true
    else
        OPT_K=21
    fi

    "$KSNP_BIN" -in "$GENOME_LIST" -outdir "$KSNP_OUT/results" -k "$OPT_K" \
        -CPU "$THREADS" 2>"$KSNP_OUT/stderr.log" || true

    # Find the resulting tree
    KSNP_TREE=""
    for candidate in "$KSNP_OUT/results/tree.ML.tre" "$KSNP_OUT/results/tree.parsimony.tre" \
                     "$KSNP_OUT/results/"*.tre; do
        if [ -f "$candidate" ]; then
            KSNP_TREE="$candidate"
            break
        fi
    done

    END=$(date +%s)
    ELAPSED=$((END - START))

    if [ -n "$KSNP_TREE" ]; then
        cp "$KSNP_TREE" "$KSNP_OUT/tree.newick"
        record_result "kSNP" "$KSNP_OUT/tree.newick" "$ELAPSED" "N/A"
    else
        echo "    [kSNP] No tree produced."
        echo -e "kSNP\tERROR\tN/A\t${ELAPSED}\tN/A\tN/A" >> "$RESULTS_TSV"
    fi
else
    echo "  [kSNP] SKIPPED (not installed)"
fi

# ============================================================
# 9. Peafowl
# ============================================================
PEAFOWL_SCRIPT="${TOOLS_DIR}/../references/Peafowl-repo/Peafowl/SingleTree/main_script.sh"
if [ -f "$PEAFOWL_SCRIPT" ] && command -v jellyfish &>/dev/null && command -v raxmlHPC-PTHREADS &>/dev/null; then
    PEAFOWL_OUT="$OUTPUT_DIR/Peafowl"
    mkdir -p "$PEAFOWL_OUT"
    echo "  [Peafowl] Running..."

    # Peafowl needs a directory with FASTA files
    PEAFOWL_IN="$PEAFOWL_OUT/input"
    mkdir -p "$PEAFOWL_IN"
    for f in "${FASTA_FILES[@]}"; do
        cp "$f" "$PEAFOWL_IN/"
    done

    # Run Peafowl
    START=$(date +%s)
    
    # Needs a dummy reference tree to run (it compares internally, but we just want the output tree)
    DUMMY_TREE="$PEAFOWL_OUT/dummy.newick"
    echo "()" > "$DUMMY_TREE"

    cd "${TOOLS_DIR}/../references/Peafowl-repo/Peafowl/SingleTree" || exit
    # Usage: ./main_script.sh <tree> <threads> <source_folder> <is_reverse_compliment>
    bash main_script.sh "$(realpath "$DUMMY_TREE")" "$THREADS" "$(realpath "$PEAFOWL_IN")" 1 > "$(realpath "$PEAFOWL_OUT")/stdout.log" 2> "$(realpath "$PEAFOWL_OUT")/stderr.log" || true
    cd "$SCRIPT_DIR/.." || exit

    # Output is saved to ${PEAFOWL_IN}_Result/Result_unrooted_tree_Kmer_*.newick
    PEAFOWL_RESULT_DIR="${PEAFOWL_IN}_Result"
    PEAFOWL_TREE=$(find "$PEAFOWL_RESULT_DIR" -name "Result_unrooted_tree_Kmer_*.newick" 2>/dev/null | head -1)

    END=$(date +%s)
    ELAPSED=$((END - START))

    if [ -n "$PEAFOWL_TREE" ] && [ -f "$PEAFOWL_TREE" ]; then
        cp "$PEAFOWL_TREE" "$PEAFOWL_OUT/tree.newick"
        record_result "Peafowl" "$PEAFOWL_OUT/tree.newick" "$ELAPSED" "N/A"
    else
        echo "    [Peafowl] No tree produced."
        echo -e "Peafowl\tERROR\tN/A\t${ELAPSED}\tN/A\tN/A" >> "$RESULTS_TSV"
    fi
    
    # Cleanup heavy temp files
    rm -rf "$PEAFOWL_IN" "$PEAFOWL_RESULT_DIR"
else
    echo "  [Peafowl] SKIPPED (not installed or missing dependencies)"
fi

echo ""
echo "============================================================"
echo "  nRF Comparison Results"
echo "============================================================"
echo ""
column -t -s $'\t' "$RESULTS_TSV"
echo ""
echo "Results saved to: $RESULTS_TSV"
echo ""
