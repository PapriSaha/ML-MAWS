/**
 * MatrixBuilder.cpp - Binary Character Matrix Construction implementation
 */

#include "MatrixBuilder.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <iomanip>

namespace mlmaws {

void MatrixBuilder::build(const std::vector<std::string>& speciesNames,
                          const std::vector<std::vector<std::string>>& mawSets,
                          bool removeConstant) {
    int m = (int)speciesNames.size();
    if (m == 0 || m != (int)mawSets.size()) {
        std::cerr << "[ERROR] Invalid input to MatrixBuilder::build" << std::endl;
        return;
    }

    speciesNames_ = speciesNames;

    // Step 1: Compute union of all MAWs (sorted, since each mawSet is sorted)
    std::set<std::string> unionSet;
    for (int i = 0; i < m; i++) {
        for (const auto& w : mawSets[i]) {
            unionSet.insert(w);
        }
    }

    std::cerr << "[INFO] Union of MAW sets: " << unionSet.size() << " unique MAWs" << std::endl;

    // Step 2: Create mapping from MAW string to column index
    std::vector<std::string> allMAWs(unionSet.begin(), unionSet.end());
    // allMAWs is already sorted (std::set iterates in sorted order)

    // Step 3: Build the binary matrix using two-pointer merge
    // matrix_[species][character] = 0 or 1
    int totalChars = (int)allMAWs.size();
    std::vector<std::vector<uint8_t>> fullMatrix(m, std::vector<uint8_t>(totalChars, 0));

    for (int i = 0; i < m; i++) {
        // Two-pointer: mawSets[i] is sorted, allMAWs is sorted
        int p = 0;
        for (int j = 0; j < totalChars && p < (int)mawSets[i].size(); j++) {
            while (p < (int)mawSets[i].size() && mawSets[i][p] < allMAWs[j]) {
                p++;
            }
            if (p < (int)mawSets[i].size() && mawSets[i][p] == allMAWs[j]) {
                fullMatrix[i][j] = 1;
                p++;
            }
        }
    }

    // Step 4: Remove constant columns (all 0 or all 1)
    if (removeConstant) {
        characterNames_.clear();
        matrix_.assign(m, std::vector<uint8_t>());

        int removedCount = 0;
        for (int j = 0; j < totalChars; j++) {
            int sum = 0;
            for (int i = 0; i < m; i++) {
                sum += fullMatrix[i][j];
            }
            // Keep only variable characters (not all 0 and not all 1)
            if (sum > 0 && sum < m) {
                characterNames_.push_back(allMAWs[j]);
                for (int i = 0; i < m; i++) {
                    matrix_[i].push_back(fullMatrix[i][j]);
                }
            } else {
                removedCount++;
            }
        }
        std::cerr << "[INFO] Removed " << removedCount
                  << " constant columns. Remaining: " << characterNames_.size() << std::endl;
    } else {
        characterNames_ = allMAWs;
        matrix_ = fullMatrix;
    }
}

void MatrixBuilder::writePhylip(const std::string& filepath) const {
    int m = numSpecies();
    int n = numCharacters();

    if (m == 0 || n == 0) {
        std::cerr << "[ERROR] Empty matrix, cannot write PHYLIP file." << std::endl;
        return;
    }

    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "[ERROR] Cannot open output file: " << filepath << std::endl;
        return;
    }

    // PHYLIP header: number of taxa and characters
    out << m << " " << n << "\n";

    // Each row: species name (padded to 10 chars or followed by spaces) then binary sequence
    for (int i = 0; i < m; i++) {
        // Truncate or pad species name to fit PHYLIP format
        std::string name = speciesNames_[i];
        // Replace problematic characters
        for (auto& c : name) {
            if (c == ' ' || c == '\t' || c == '(' || c == ')' ||
                c == '[' || c == ']' || c == ':' || c == ';' || c == ',') {
                c = '_';
            }
        }
        // Truncate to 50 chars max for RAxML compatibility
        if (name.size() > 50) name = name.substr(0, 50);

        out << name << "  ";
        for (int j = 0; j < n; j++) {
            out << (int)matrix_[i][j];
        }
        out << "\n";
    }

    out.close();
    std::cerr << "[INFO] PHYLIP matrix written: " << filepath
              << " (" << m << " taxa, " << n << " characters)" << std::endl;
}

void MatrixBuilder::writeTSV(const std::string& filepath) const {
    int m = numSpecies();
    int n = numCharacters();

    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "[ERROR] Cannot open output file: " << filepath << std::endl;
        return;
    }

    // Header row: species names
    for (int i = 0; i < m; i++) {
        if (i > 0) out << "\t";
        out << speciesNames_[i];
    }
    out << "\n";

    // Data rows: one per MAW character (transposed for entropy computation)
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            if (i > 0) out << "\t";
            out << (int)matrix_[i][j];
        }
        out << "\n";
    }

    out.close();
}

void MatrixBuilder::capCharacters(int maxChars) {
    int m = numSpecies();
    int n = numCharacters();

    if (n <= maxChars || maxChars <= 0) return;

    // Score each character by parsimony-informativeness:
    // Best = closest to 50% presence (highest variance)
    // Score = min(sum, m - sum) → higher is more informative
    std::vector<std::pair<int, int>> scores(n); // (score, original_index)
    for (int j = 0; j < n; j++) {
        int sum = 0;
        for (int i = 0; i < m; i++) {
            sum += matrix_[i][j];
        }
        int informativeness = std::min(sum, m - sum);
        scores[j] = {informativeness, j};
    }

    // Sort by informativeness (descending), then by index (ascending) for stability
    std::sort(scores.begin(), scores.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });

    // Keep top maxChars, then re-sort by original index to preserve column order
    scores.resize(maxChars);
    std::sort(scores.begin(), scores.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    // Rebuild matrix and character names
    std::vector<std::string> newCharNames;
    std::vector<std::vector<uint8_t>> newMatrix(m);

    for (const auto& [score, idx] : scores) {
        newCharNames.push_back(characterNames_[idx]);
        for (int i = 0; i < m; i++) {
            newMatrix[i].push_back(matrix_[i][idx]);
        }
    }

    characterNames_ = std::move(newCharNames);
    matrix_ = std::move(newMatrix);

    std::cerr << "[INFO] Capped to " << numCharacters() << " most informative characters" << std::endl;
}

} // namespace mlmaws
