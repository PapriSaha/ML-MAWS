#!/bin/bash
# ============================================================
# install_baselines.sh - Install all baseline phylogeny methods
#
# Installs 8 alignment-free methods for nRF comparison:
#   1. Mash (MinHash distances)
#   2. Skmer (k-mer sketching)
#   3. CAFE (k-mer composition)
#   4. FFP (Feature Frequency Profile)
#   5. Co-phylog (composition distance)
#   6. andi (anchor distance)
#   7. CD-MAWS-SA (composition distance on MAWs)
#   8. Peafowl (ML on k-mer matrix)
#
# !bash scripts/install_baselines.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}/.."
TOOLS_DIR="${PROJECT_DIR}/tools"
mkdir -p "$TOOLS_DIR"

echo ""
echo "============================================================"
echo "  Installing Baseline Methods"
echo "============================================================"
echo ""

# ---- 1. Mash ----
echo "[1/8] Mash..."
if command -v mash &>/dev/null; then
    echo "  Already installed: $(mash --version 2>&1 | head -1)"
else
    conda install -y -q -c bioconda mash 2>/dev/null || {
        cd "$TOOLS_DIR"
        wget -q https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
        tar xf mash-Linux64-v2.3.tar
        cp mash-Linux64-v2.3/mash /usr/local/bin/
        rm -rf mash-Linux64-v2.3*
        echo "  Installed Mash."
    }
fi

# ---- 2. Skmer ----
echo "[2/8] Skmer..."
if command -v skmer &>/dev/null; then
    echo "  Already installed."
else
    pip install -q skmer 2>/dev/null || {
        cd "$TOOLS_DIR"
        git clone --depth 1 https://github.com/shahab-sarmashghi/Skmer.git 2>/dev/null || true
        cd Skmer && pip install -q . 2>/dev/null || true
        echo "  Installed Skmer."
    }
fi

# ---- 3. CAFE ----
echo "[3/8] CAFE..."
if [ -f "$TOOLS_DIR/CAFE/cafe" ] || command -v cafe &>/dev/null; then
    echo "  Already installed."
else
    cd "$TOOLS_DIR"
    git clone --depth 1 https://github.com/LuoGroup2023/CAFE.git 2>/dev/null || true
    if [ -d "CAFE" ]; then
        cd CAFE
        make 2>/dev/null || g++ -O2 -o cafe *.cpp 2>/dev/null || true
        echo "  Installed CAFE."
    fi
fi

# ---- 4. FFP (Feature Frequency Profile) ----
echo "[4/8] FFP..."
if command -v ffpry &>/dev/null || [ -f "$TOOLS_DIR/ffp/ffpry" ]; then
    echo "  Already installed."
else
    conda install -y -q -c bioconda ffp 2>/dev/null || {
        cd "$TOOLS_DIR"
        git clone --depth 1 https://github.com/apetkau/ffp-3.19-custom.git ffp 2>/dev/null || true
        if [ -d "ffp" ]; then
            cd ffp && ./configure --disable-gui && make 2>/dev/null || true
            echo "  Installed FFP."
        fi
    }
fi

# ---- 5. Co-phylog ----
echo "[5/8] Co-phylog..."
if [ -f "$TOOLS_DIR/co-phylog/co-phylog" ] || command -v co-phylog &>/dev/null; then
    echo "  Already installed."
else
    cd "$TOOLS_DIR"
    git clone --depth 1 https://github.com/coPhylogDev/co-phylog.git 2>/dev/null || true
    if [ -d "co-phylog" ]; then
        cd co-phylog
        make 2>/dev/null || true
        echo "  Installed Co-phylog."
    fi
fi

# ---- 5b. kSNP4 (replaces kSNP3) ----
echo "[5b/8] kSNP4..."
if command -v kSNP4 &>/dev/null || [ -f "$TOOLS_DIR/kSNP4/kSNP4" ]; then
    echo "  Already installed."
else
    cd "$TOOLS_DIR"
    curl -L "https://sourceforge.net/projects/ksnp/files/kSNP4.1pkg/kSNP4.1_Linux_package.zip/download" -o kSNP4.zip 2>/dev/null || true
    if [ -f "kSNP4.zip" ]; then
        unzip -qo kSNP4.zip -d kSNP4 2>/dev/null || true
        # Find the actual binary directory
        KSNP_BIN=$(find kSNP4 -name "kSNP4" -type f 2>/dev/null | head -1)
        if [ -n "$KSNP_BIN" ]; then
            KSNP_DIR=$(dirname "$KSNP_BIN")
            export PATH="$PATH:$KSNP_DIR"
            echo "export PATH=\$PATH:$KSNP_DIR" >> ~/.bashrc
            echo "  Installed kSNP4."
        else
            echo "  [WARN] kSNP4 binary not found after extraction."
        fi
        rm -f kSNP4.zip
    else
        echo "  [WARN] kSNP4 download failed."
    fi
fi

# ---- 6. andi ----
echo "[6/8] andi..."
if command -v andi &>/dev/null; then
    echo "  Already installed."
else
    conda install -y -q -c bioconda andi 2>/dev/null || {
        apt-get install -y -qq libdivsufsort-dev libgsl-dev 2>/dev/null || true
        cd "$TOOLS_DIR"
        git clone --depth 1 https://github.com/EvolBioInf/andi.git 2>/dev/null || true
        if [ -d "andi" ]; then
            cd andi
            autoreconf -fi 2>/dev/null && ./configure && make 2>/dev/null || true
            echo "  Installed andi."
        fi
    }
fi

# ---- 7. CD-MAWS-SA ----
echo "[7/8] CD-MAWS-SA..."
CDMAWS_REPO="${TOOLS_DIR}/cd-maws-sa"
CDMAWS_DIR="${CDMAWS_REPO}/cd-maws-sa"
if [ ! -d "$CDMAWS_REPO" ]; then
    cd "$TOOLS_DIR"
    # Try public clone first; fallback to Drive copy
    git clone --depth 1 https://github.com/TamimEhsan/cd-maws-sa.git cd-maws-sa 2>/dev/null || \
    cp -r /content/drive/MyDrive/ML-MAWS-Data/cd-maws-sa "$CDMAWS_REPO" 2>/dev/null || true
fi
if [ -d "$CDMAWS_DIR" ] && [ -f "$CDMAWS_DIR/maw-gen.cpp" ]; then
    cd "$CDMAWS_DIR"
    g++ -O2 -std=c++17 -o maw-gen maw-gen.cpp 2>/dev/null && echo "  Compiled maw-gen." || echo "  [WARN] maw-gen compile failed."
    g++ -O2 -std=c++17 -o dist-calc dist-calc.cpp 2>/dev/null && echo "  Compiled dist-calc." || echo "  [WARN] dist-calc compile failed."
else
    echo "  [WARN] cd-maws-sa source not found."
    echo "         Options:"
    echo "           1. Make https://github.com/TamimEhsan/cd-maws-sa public"
    echo "           2. Upload source to Drive at: ML-MAWS-Data/cd-maws-sa/"
fi

# ---- 8. Peafowl (requires Jellyfish and RAxML) ----
echo "[8/8] Peafowl dependencies..."
if command -v jellyfish &>/dev/null; then
    echo "  Jellyfish already installed."
else
    conda install -y -q -c bioconda jellyfish 2>/dev/null || {
        apt-get install -y -qq jellyfish 2>/dev/null || true
    }
fi

if command -v raxmlHPC-PTHREADS &>/dev/null; then
    echo "  RAxML already installed."
else
    conda install -y -q -c bioconda raxml 2>/dev/null || {
        apt-get install -y -qq raxml 2>/dev/null || true
    }
fi
PEAFOWL_DIR="${PROJECT_DIR}/references/Peafowl-repo"
if [ ! -d "$PEAFOWL_DIR" ]; then
    mkdir -p "${PROJECT_DIR}/references"
    cd "${PROJECT_DIR}/references"
    git clone https://github.com/hasin-abrar/Peafowl-repo.git 2>/dev/null || true
fi

if [ -d "$PEAFOWL_DIR/Peafowl/SingleTree" ]; then
    cd "$PEAFOWL_DIR/Peafowl/SingleTree"
    g++ -O2 -o kmerMerge kmerMerge.cpp -lpthread 2>/dev/null || true
    g++ -O2 -o entropy entropy.cpp 2>/dev/null || true
    echo "  Compiled Peafowl components."
else
    echo "  [WARN] Peafowl source not found after clone attempt."
fi

echo ""
echo "============================================================"
echo "  Verification"
echo "============================================================"
echo ""
echo -n "  mash:       "; command -v mash >/dev/null && echo "OK" || echo "MISSING"
echo -n "  skmer:      "; command -v skmer >/dev/null && echo "OK" || echo "MISSING"
echo -n "  andi:       "; command -v andi >/dev/null && echo "OK" || echo "MISSING"
echo -n "  jellyfish:  "; command -v jellyfish >/dev/null && echo "OK" || echo "MISSING"
echo -n "  cd-maws-sa: "; [ -f "$CDMAWS_DIR/maw-gen" ] && echo "OK" || echo "MISSING"
echo ""
echo "Done."
echo ""
