#!/bin/bash
# ML-MAWS
# Maximum Likelihood Phylogeny from Minimal Absent Word Characters
# Prerequisites:
#   - ml-maws binary (compile with: make)
#   - RAxML (raxmlHPC or raxmlHPC-PTHREADS) OR IQ-TREE (iqtree2)
#   - Python 3 with dendropy

set -e

# Default parameters
INPUT_DIR=""
OUTPUT_DIR="ml_maws_results"
THREADS=4
STRAND_FLAG=""
ML_TOOL="raxml"    # raxml or iqtree
FIXED_LENGTH=0      # 0 = auto-select
LMIN=3
LMAX=12
REFERENCE_TREE=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BINARY="${SCRIPT_DIR}/../ml-maws"

# ---- Parse arguments ----
usage() {
    echo ""
    echo "ML-MAWS Pipeline Script"
    echo "======================="
    echo ""
    echo "Usage: $0 -i <input_dir> [options]"
    echo ""
    echo "Required:"
    echo "  -i <path>       Input directory of FASTA files"
    echo ""
    echo "Optional:"
    echo "  -o <path>       Output directory (default: ml_maws_results)"
    echo "  -t <int>        Number of threads (default: 4)"
    echo "  -l <int>        Fixed MAW length (0 = auto, default: 0)"
    echo "  --lmin <int>    Min MAW length for search (default: 3)"
    echo "  --lmax <int>    Max MAW length for search (default: 12)"
    echo "  --strand        Enable strand-aware MAW filtering"
    echo "  --iqtree        Use IQ-TREE instead of RAxML"
    echo "  --ref <tree>    Reference tree for nRF comparison"
    echo "  -h              Show this help"
    echo ""
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)       INPUT_DIR="$2"; shift 2 ;;
        -o)       OUTPUT_DIR="$2"; shift 2 ;;
        -t)       THREADS="$2"; shift 2 ;;
        -l)       FIXED_LENGTH="$2"; shift 2 ;;
        --lmin)   LMIN="$2"; shift 2 ;;
        --lmax)   LMAX="$2"; shift 2 ;;
        --strand) STRAND_FLAG="--strand"; shift ;;
        --iqtree) ML_TOOL="iqtree"; shift ;;
        --ref)    REFERENCE_TREE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *)        echo "Unknown option: $1"; usage ;;
    esac
done

if [ -z "$INPUT_DIR" ]; then
    echo "Error: No input directory specified."
    usage
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

# ---- Check prerequisites ----
if [ ! -f "$BINARY" ]; then
    echo "ml-maws binary not found. Compiling..."
    cd "${SCRIPT_DIR}/.."
    make
    cd -
fi

echo ""
echo ""
echo "  Input:    $INPUT_DIR"
echo "  Output:   $OUTPUT_DIR"
echo "  Threads:  $THREADS"
echo "  ML Tool:  $ML_TOOL"
echo "  Strand:   ${STRAND_FLAG:-disabled}"
echo "  Length:    ${FIXED_LENGTH} (0=auto)"
echo ""

# ---- Build ml-maws arguments ----
ARGS="-i $INPUT_DIR -o $OUTPUT_DIR --threads $THREADS --lmin $LMIN --lmax $LMAX"

if [ "$FIXED_LENGTH" -ne 0 ]; then
    ARGS="$ARGS -l $FIXED_LENGTH"
fi

if [ -n "$STRAND_FLAG" ]; then
    ARGS="$ARGS --strand"
fi

if [ "$ML_TOOL" == "iqtree" ]; then
    ARGS="$ARGS --iqtree"
fi

# ---- Run ML-MAWS ----
echo "Running: $BINARY $ARGS"
echo ""

START_TIME=$(date +%s)
$BINARY $ARGS
END_TIME=$(date +%s)

ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "ML-MAWS completed in ${ELAPSED} seconds."

# ---- Compare with reference tree (if provided) ----
if [ -n "$REFERENCE_TREE" ] && [ -f "$REFERENCE_TREE" ]; then
    RESULT_TREE="$OUTPUT_DIR/ML_MAWS_tree.newick"
    if [ -f "$RESULT_TREE" ]; then
        echo ""
        echo "Comparing with reference tree..."
        python3 "${SCRIPT_DIR}/compare_trees.py" "$REFERENCE_TREE" "$RESULT_TREE"
    else
        echo "Warning: Result tree not found at $RESULT_TREE"
    fi
fi

echo ""
echo "Done. Results in: $OUTPUT_DIR/"
echo ""
