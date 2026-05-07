/**
 * MatrixBuilder.h - Binary Character Matrix Construction
 * 
 * Builds the m × |U| binary matrix where:
 *   M[i][w] = 1 if MAW w is present in species i's MAW set
 *   M[i][w] = 0 otherwise
 * 
 * Outputs in PHYLIP format for RAxML / IQ-TREE.
 */

#ifndef MATRIX_BUILDER_H
#define MATRIX_BUILDER_H

#include <string>
#include <vector>
#include <map>
#include <set>

namespace mlmaws {

/**
 * Builds and manages the binary MAW character matrix.
 */
class MatrixBuilder {
public:
    MatrixBuilder() = default;
    ~MatrixBuilder() = default;

    /**
     * Build the binary character matrix from per-species MAW sets.
     * 
     * @param speciesNames  Names of each species
     * @param mawSets       MAW sets for each species (sorted vectors)
     * @param removeConstant If true, remove columns that are all-0 or all-1
     */
    void build(const std::vector<std::string>& speciesNames,
               const std::vector<std::vector<std::string>>& mawSets,
               bool removeConstant = true);

    /**
     * Write the matrix in PHYLIP format (binary characters).
     * Suitable for RAxML (-m BINGAMMA) or IQ-TREE (-st BIN).
     */
    void writePhylip(const std::string& filepath) const;

    /**
     * Write the matrix in tab-separated format (for entropy calculation).
     * First row: species names. Each subsequent row: one MAW character.
     */
    void writeTSV(const std::string& filepath) const;

    /** Get number of species (rows). */
    int numSpecies() const { return (int)speciesNames_.size(); }

    /** Get number of characters (columns) after filtering. */
    int numCharacters() const { return (int)characterNames_.size(); }

    /** Get the full binary matrix. */
    const std::vector<std::vector<uint8_t>>& getMatrix() const { return matrix_; }

    /** Get species names. */
    const std::vector<std::string>& getSpeciesNames() const { return speciesNames_; }

    /** Get character (MAW) names. */
    const std::vector<std::string>& getCharacterNames() const { return characterNames_; }

    /**
     * Cap the number of characters by keeping the most parsimony-informative ones.
     * Characters are ranked by variance (distance from 0.5 frequency is worst).
     * @param maxChars  Maximum number of characters to keep
     */
    void capCharacters(int maxChars);

private:
    std::vector<std::string> speciesNames_;
    std::vector<std::string> characterNames_;  // the MAW strings (union)
    std::vector<std::vector<uint8_t>> matrix_; // [species][character]
};

} // namespace mlmaws

#endif // MATRIX_BUILDER_H
