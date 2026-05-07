#!/bin/bash
# run_all_benchmarks.sh 

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}/.."
ML_MAWS_BIN="${PROJECT_DIR}/ml-maws"
RESULTS_DIR="${PROJECT_DIR}/benchmark/results"
FIGURES_DIR="${PROJECT_DIR}/benchmark/figures"
THREADS=2
MAX_CHARS=50000

DRIVE_DATA="/content/drive/MyDrive/ML-MAWS-Data"

echo ""
echo "============================================================"
echo "  ML-MAWS Benchmark"
echo "============================================================"
echo ""

# ---- Build ----
echo "[1] Building ML-MAWS..."
cd "$PROJECT_DIR"
if [ ! -f ml-maws ]; then
    make clean && make
fi
echo "  Build OK"
echo ""

# ---- Setup ----
mkdir -p "$RESULTS_DIR" "$FIGURES_DIR"

TIMING_JSON="${RESULTS_DIR}/timing_results.json"
echo '{}' > "$TIMING_JSON"

# ---- Helper: save timing to JSON ----
save_timing() {
    local dataset="$1"
    local method="$2"
    local elapsed="$3"

    python3 -c "
import json, os
fpath = '$TIMING_JSON'
data = {}
if os.path.exists(fpath):
    with open(fpath) as f:
        data = json.load(f)
if '$dataset' not in data:
    data['$dataset'] = {}
data['$dataset']['$method'] = {
    'time_seconds': $elapsed
}
with open(fpath, 'w') as f:
    json.dump(data, f, indent=2)
" 2>/dev/null || true
}

# ============================================================
# Define datasets: name -> Drive path
# ============================================================
# Directory-based datasets (one .fasta per species)
declare -A DATASETS
DATASETS["fish-mtdna"]="$DRIVE_DATA/fish_mtDNA"
DATASETS["ecoli-shigella"]="$DRIVE_DATA/ecoli-shigella"
DATASETS["ecoli-shigella-hgt"]="$DRIVE_DATA/ecoli_shigella(HGT)"
DATASETS["yersinia-hgt"]="$DRIVE_DATA/yersinia-HGT"
# Plants skipped: 4.8GB total, exceeds Colab RAM
# DATASETS["plants"]="$DRIVE_DATA/Plants"

# Simulated HGT has sub-datasets
declare -A SIM_DATASETS
SIM_DATASETS["simulated-hgt-0"]="$DRIVE_DATA/simulated-HGT/hgt_0"
SIM_DATASETS["simulated-hgt-250"]="$DRIVE_DATA/simulated-HGT/hgt_250"
SIM_DATASETS["simulated-hgt-500"]="$DRIVE_DATA/simulated-HGT/hgt_500"
SIM_DATASETS["simulated-hgt-750"]="$DRIVE_DATA/simulated-HGT/hgt_750"
SIM_DATASETS["simulated-hgt-1000"]="$DRIVE_DATA/simulated-HGT/hgt_1000"

# Single multi-FASTA files (all species in one file, each >header = one species)
declare -A FILE_DATASETS
FILE_DATASETS["coronavirus"]="$DRIVE_DATA/coronavirus.fasta"
FILE_DATASETS["ebolavirus"]="$DRIVE_DATA/ebolavirus.fasta"
FILE_DATASETS["influenza"]="$DRIVE_DATA/influenza.fasta"
FILE_DATASETS["mammal-mtdna"]="$DRIVE_DATA/mammal_mtdna.fasta"
FILE_DATASETS["rhinovirus"]="$DRIVE_DATA/rhinovirus.fasta"

echo "[2] Checking datasets..."
echo ""

# ---- Check which datasets exist ----
VALID_DATASETS=()
for ds in "${!DATASETS[@]}"; do
    path="${DATASETS[$ds]}"
    count=$(ls "$path"/*.fasta 2>/dev/null | wc -l)
    if [ "$count" -gt 0 ]; then
        echo "  ✓ $ds: $count species ($path)"
        VALID_DATASETS+=("$ds")
    else
        echo "  ✗ $ds: NOT FOUND ($path)"
    fi
done

for ds in "${!SIM_DATASETS[@]}"; do
    path="${SIM_DATASETS[$ds]}"
    count=$(ls "$path"/*.fasta 2>/dev/null | wc -l)
    if [ "$count" -gt 0 ]; then
        echo "  ✓ $ds: $count species ($path)"
        VALID_DATASETS+=("$ds")
    else
        echo "  ✗ $ds: NOT FOUND ($path)"
    fi
done

for ds in "${!FILE_DATASETS[@]}"; do
    path="${FILE_DATASETS[$ds]}"
    if [ -f "$path" ] && [ -s "$path" ]; then
        species=$(grep -c "^>" "$path" 2>/dev/null || echo "0")
        echo "  ✓ $ds: $species species ($path)"
        VALID_DATASETS+=("$ds")
    else
        echo "  ✗ $ds: NOT FOUND ($path)"
    fi
done

if [ ${#VALID_DATASETS[@]} -eq 0 ]; then
    echo ""
    echo "[ERROR] No datasets found! Check that Drive is mounted and paths are correct."
    exit 1
fi

echo ""
echo "[3] Running ML-MAWS on ${#VALID_DATASETS[@]} datasets..."
echo ""

# ---- Run on each dataset ----
for ds in "${VALID_DATASETS[@]}"; do
    # Get input path
    if [ -n "${DATASETS[$ds]+x}" ]; then
        INPUT_DIR="${DATASETS[$ds]}"
        SPECIES_COUNT=$(ls "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l)
    elif [ -n "${SIM_DATASETS[$ds]+x}" ]; then
        INPUT_DIR="${SIM_DATASETS[$ds]}"
        SPECIES_COUNT=$(ls "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l)
    else
        INPUT_DIR="${FILE_DATASETS[$ds]}"
        SPECIES_COUNT=$(grep -c "^>" "$INPUT_DIR" 2>/dev/null || echo "?")
    fi

    echo "============================================"
    echo "  Dataset: $ds"
    echo "  Input:   $INPUT_DIR"
    echo "  Species: $SPECIES_COUNT"
    echo "============================================"

    # --- Variant 1: ML-MAWS + strand + IQ-TREE ---
    OUT="$RESULTS_DIR/$ds/ml-maws-iqtree"
    mkdir -p "$OUT"
    echo "  Running: ML-MAWS-IQTree..."
    SECONDS=0
    "$ML_MAWS_BIN" -i "$INPUT_DIR" -o "$OUT" --strand --iqtree --max-chars $MAX_CHARS --threads $THREADS 2>"$OUT/stderr.log" || true
    ELAPSED=$SECONDS
    echo "    Time: ${ELAPSED}s"
    save_timing "$ds" "ML-MAWS-IQTree" "$ELAPSED"

    # --- Variant 2: ML-MAWS + strand + RAxML ---
    OUT2="$RESULTS_DIR/$ds/ml-maws-raxml"
    mkdir -p "$OUT2"
    echo "  Running: ML-MAWS-RAxML..."
    SECONDS=0
    "$ML_MAWS_BIN" -i "$INPUT_DIR" -o "$OUT2" --strand --max-chars $MAX_CHARS --threads $THREADS 2>"$OUT2/stderr.log" || true
    ELAPSED=$SECONDS
    echo "    Time: ${ELAPSED}s"
    save_timing "$ds" "ML-MAWS-RAxML" "$ELAPSED"

    # --- Variant 3: ML-MAWS no-strand + IQ-TREE (ablation) ---
    OUT3="$RESULTS_DIR/$ds/ml-maws-nostrand"
    mkdir -p "$OUT3"
    echo "  Running: ML-MAWS-NoStrand..."
    SECONDS=0
    "$ML_MAWS_BIN" -i "$INPUT_DIR" -o "$OUT3" --iqtree --max-chars $MAX_CHARS --threads $THREADS 2>"$OUT3/stderr.log" || true
    ELAPSED=$SECONDS
    echo "    Time: ${ELAPSED}s"
    save_timing "$ds" "ML-MAWS-NoStrand" "$ELAPSED"

    # Copy reference tree if it exists (for nRF evaluation)
    DS_RESULT="$RESULTS_DIR/$ds"
    REF_FOUND=false

    # Search in multiple locations and patterns
    for ref_candidate in "$INPUT_DIR/../reference.newick" "$INPUT_DIR/../ref.newick" \
                         "$INPUT_DIR/../true_tree.newick" "$INPUT_DIR/reference.newick"; do
        if [ -f "$ref_candidate" ]; then
            cp "$ref_candidate" "$DS_RESULT/reference.newick"
            echo "  Reference tree: copied from $(realpath "$ref_candidate")"
            REF_FOUND=true
            break
        fi
    done

    # Also check inside a reference/ subfolder (any .newick file)
    if [ "$REF_FOUND" = false ]; then
        for ref_dir in "$INPUT_DIR/reference" "$INPUT_DIR/../reference" \
                       "$INPUT_DIR/ref" "$INPUT_DIR/../ref"; do
            if [ -d "$ref_dir" ]; then
                REF_FILE=$(find "$ref_dir" -maxdepth 1 -name "*.newick" -o -name "*.nwk" -o -name "*.tree" | head -1)
                if [ -n "$REF_FILE" ]; then
                    cp "$REF_FILE" "$DS_RESULT/reference.newick"
                    echo "  Reference tree: copied from $REF_FILE"
                    REF_FOUND=true
                    break
                fi
            fi
        done
    fi

    if [ "$REF_FOUND" = false ]; then
        echo "  Reference tree: not found (nRF evaluation will be skipped)"
    fi

    echo ""
done

# ---- Collect entropy curves ----
echo "[4] Collecting results..."
ENTROPY_DIR="${RESULTS_DIR}/entropy_curves"
mkdir -p "$ENTROPY_DIR"
for ds in "${VALID_DATASETS[@]}"; do
    for variant in ml-maws-iqtree ml-maws-raxml ml-maws-nostrand; do
        ENTROPY_FILE="$RESULTS_DIR/$ds/$variant/entropy_results.tsv"
        if [ -f "$ENTROPY_FILE" ]; then
            cp "$ENTROPY_FILE" "$ENTROPY_DIR/${ds}_${variant}_entropy.tsv"
        fi
    done
done

# ---- Evaluate nRF/nQD ----
EVAL_JSON="${RESULTS_DIR}/evaluation_results.json"
echo ""
echo "[5] Computing nRF and nQD metrics..."
python3 "${SCRIPT_DIR}/evaluate.py" --results-dir "$RESULTS_DIR" --output "$EVAL_JSON" 2>/dev/null || true
echo ""


echo "[6] Generating figures..."

# nRF and nQD comparison
python3 "${SCRIPT_DIR}/plot_results.py" --results "$EVAL_JSON" --output-dir "$FIGURES_DIR" 2>/dev/null || true

# Timing and memory plots
python3 "${SCRIPT_DIR}/plot_results.py" --timing "$TIMING_JSON" --output-dir "$FIGURES_DIR" 2>/dev/null || true

# Entropy curves
python3 "${SCRIPT_DIR}/plot_results.py" --entropy-dir "$ENTROPY_DIR" --output-dir "$FIGURES_DIR" 2>/dev/null || true

# Tree visualizations
for ds in "${VALID_DATASETS[@]}"; do
    TREE="$RESULTS_DIR/$ds/ml-maws-iqtree/ML_MAWS_tree.newick"
    if [ -f "$TREE" ]; then
        python3 "${SCRIPT_DIR}/plot_results.py" --tree "$TREE" \
            --tree-title "ML-MAWS: $ds" --output-dir "$FIGURES_DIR" 2>/dev/null || true
        mv "$FIGURES_DIR/tree.pdf" "$FIGURES_DIR/tree_${ds}.pdf" 2>/dev/null || true
    fi
done

# Complexity breakdown
python3 "${SCRIPT_DIR}/plot_results.py" --complexity-dir "$RESULTS_DIR" --output-dir "$FIGURES_DIR" 2>/dev/null || true

echo ""
echo "============================================================"
echo "  Benchmark Complete!"
echo "============================================================"
echo ""
echo "  Results:  $RESULTS_DIR/"
echo "  Timing:   $TIMING_JSON"
echo "  Figures:  $FIGURES_DIR/"
echo ""
echo "  Trees generated:"
for ds in "${VALID_DATASETS[@]}"; do
    TREE="$RESULTS_DIR/$ds/ml-maws-iqtree/ML_MAWS_tree.newick"
    if [ -f "$TREE" ]; then
        echo "    ✓ $ds: $TREE"
    else
        echo "    ✗ $ds: no tree (check $RESULTS_DIR/$ds/ml-maws-iqtree/stderr.log)"
    fi
done
echo ""
echo "  Copy results to Drive:"
echo "    cp -r $RESULTS_DIR /content/drive/MyDrive/ML-MAWS/benchmark/"
echo "    cp -r $FIGURES_DIR /content/drive/MyDrive/ML-MAWS/benchmark/"
echo ""
