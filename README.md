# ML-MAWS

A phylogenetic tree construction method using Alignment-free maximum likelihood on Minimal Absent Words.

## Build

```bash
git clone https://github.com/PapriSaha/ML-MAWS.git
cd ML-MAWS
make
```

Requires GCC 11+ (C++17, OpenMP).

## Usage

```bash
# Basic (auto-selects optimal MAW length)
./ml-maws -i <input_fasta_dir> -o <output_dir>

# With strand-aware filtering + IQ-TREE
./ml-maws -i data/ecoli -o results/ecoli --strand --iqtree

# Fixed MAW length, skip tree estimation
./ml-maws -i data/fish -o results/fish -l 8 --no-raxml
```

## Benchmarking

```bash
# Install baseline methods
bash scripts/install_baselines.sh

# Run ML-MAWS + baselines on one dataset
bash scripts/run_baselines.sh \
  -i <input_dir> -o <output_dir> \
  --ref <reference.newick> -t 4

# Run all datasets
bash scripts/run_all_benchmarks.sh
```

## Dependencies

- **GCC 11+** (C++17, OpenMP)
- [RAxML](https://github.com/stamatak/standard-RAxML) or [IQ-TREE](http://www.iqtree.org/)
