/**
 * EntropySelector.cpp - Entropy-based optimal MAW length selection implementation
 */

#include "EntropySelector.h"
#include <iostream>
#include <algorithm>

namespace mlmaws {

double EntropySelector::computeEntropy(const std::vector<std::vector<uint8_t>>& matrix, int m) {
    if (matrix.empty() || matrix[0].empty()) return 0.0;

    int n = (int)matrix[0].size();
    double totalEntropy = 0.0;

    for (int j = 0; j < n; j++) {
        int ones = 0;
        for (int i = 0; i < m; i++) {
            ones += matrix[i][j];
        }
        if (ones == 0 || ones == m) continue; // constant column

        double p1 = (double)ones / m;
        double p0 = 1.0 - p1;
        totalEntropy += -(p0 * std::log2(p0) + p1 * std::log2(p1));
    }

    return totalEntropy;
}

int EntropySelector::selectBestLength(const std::vector<EntropyResult>& results) {
    if (results.empty()) return -1;

    // First pass: find the actual best
    int bestLength = results[0].length;
    double bestEntropy = results[0].entropy;
    for (const auto& r : results) {
        if (r.entropy > bestEntropy) {
            bestEntropy = r.entropy;
            bestLength = r.length;
        }
    }

    // Second pass: print table with correct marker
    std::cerr << "\n[INFO] === Entropy-based Length Selection ===" << std::endl;
    std::cerr << "  Length | Characters | Entropy" << std::endl;
    std::cerr << "  -------|------------|--------" << std::endl;

    for (const auto& r : results) {
        std::cerr << "     " << r.length
                  << "  |   " << r.numCharacters
                  << "    | " << r.entropy;
        if (r.length == bestLength) {
            std::cerr << "  <-- best";
        }
        std::cerr << std::endl;
    }

    std::cerr << "[INFO] Selected optimal length: " << bestLength
              << " (entropy = " << bestEntropy << ")" << std::endl;

    return bestLength;
}

std::vector<int> EntropySelector::selectTopLengths(
    const std::vector<EntropyResult>& results, int topK, int minChars) {

    if (results.empty()) return {};

    // Sort by entropy (descending), keep only those with enough characters
    std::vector<EntropyResult> candidates;
    for (const auto& r : results) {
        if (r.numCharacters >= minChars && r.entropy > 0.0) {
            candidates.push_back(r);
        }
    }

    // If too few candidates, relax the minChars requirement
    if ((int)candidates.size() < topK) {
        candidates.clear();
        for (const auto& r : results) {
            if (r.numCharacters > 0 && r.entropy > 0.0) {
                candidates.push_back(r);
            }
        }
    }

    // Sort by entropy descending
    std::sort(candidates.begin(), candidates.end(),
              [](const EntropyResult& a, const EntropyResult& b) {
                  return a.entropy > b.entropy;
              });

    // Take topK
    std::vector<int> selected;
    for (int i = 0; i < std::min(topK, (int)candidates.size()); i++) {
        selected.push_back(candidates[i].length);
    }
    std::sort(selected.begin(), selected.end());

    std::cerr << "[INFO] Multi-length aggregation: selected lengths = {";
    for (size_t i = 0; i < selected.size(); i++) {
        if (i > 0) std::cerr << ", ";
        std::cerr << selected[i];
    }
    std::cerr << "}" << std::endl;

    return selected;
}

void EntropySelector::computeAdaptiveRange(int avgSeqLen, int numSpecies,
                                            int& lmin, int& lmax) {
    // Adaptive strategy based on input characteristics
    //
    // Key insight: MAW count grows with sequence length.
    // Short sequences need shorter MAWs to have enough signal.
    // Large species counts need more discriminative (longer) MAWs.

    if (avgSeqLen < 500) {
        // Very short sequences: viral segments, gene sequences
        lmin = 2;
        lmax = 6;
    } else if (avgSeqLen < 2000) {
        // Short sequences: mitochondrial genes, small genomes
        lmin = 2;
        lmax = 8;
    } else if (avgSeqLen < 50000) {
        // Medium sequences: mitochondrial genomes, plasmids
        lmin = 3;
        lmax = 10;
    } else if (avgSeqLen < 500000) {
        // Long sequences: bacterial genomes
        lmin = 4;
        lmax = 12;
    } else {
        // Very long: eukaryotic chromosomes
        lmin = 5;
        lmax = 14;
    }

    // For large species counts, expand range for better discrimination
    if (numSpecies > 50) {
        lmax = std::min(lmax + 2, 16);
    }
    if (numSpecies > 100) {
        lmax = std::min(lmax + 2, 18);
    }

    std::cerr << "[INFO] Adaptive length range: [" << lmin << ", " << lmax << "]"
              << " (avgSeqLen=" << avgSeqLen << ", nSpecies=" << numSpecies << ")"
              << std::endl;
}

} // namespace mlmaws
