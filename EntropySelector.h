/**
 * EntropySelector.h - Entropy-based optimal MAW length selection
 * 
 * Inspired by Peafowl's k-entropy approach (Zahin et al., BMC Bioinf. 2025).
 * Selects the MAW length that maximizes cumulative Shannon entropy across
 * the binary character matrix columns.
 */

#ifndef ENTROPY_SELECTOR_H
#define ENTROPY_SELECTOR_H

#include <vector>
#include <string>
#include <cmath>

namespace mlmaws {

/**
 * Result of entropy computation for a single MAW length.
 */
struct EntropyResult {
    int length;           // MAW length tested
    double entropy;       // cumulative Shannon entropy
    int numCharacters;    // number of variable characters at this length
};

/**
 * Selects the optimal MAW length based on Shannon entropy of the
 * resulting binary character matrix.
 */
class EntropySelector {
public:
    EntropySelector() = default;
    ~EntropySelector() = default;

    /**
     * Compute entropy for a binary matrix (given as column-wise presence counts).
     * 
     * @param matrix   Binary matrix [species][character]
     * @param m        Number of species
     * @return         Cumulative Shannon entropy
     */
    static double computeEntropy(const std::vector<std::vector<uint8_t>>& matrix, int m);

    /**
     * Compute entropy for a single column.
     */
    static double columnEntropy(const std::vector<uint8_t>& column, int m) {
        int ones = 0;
        for (auto v : column) ones += v;
        if (ones == 0 || ones == m) return 0.0;
        double p1 = (double)ones / m;
        double p0 = 1.0 - p1;
        return -(p0 * std::log2(p0) + p1 * std::log2(p1));
    }

    /**
     * Select the optimal MAW length from a range.
     * @param results  Entropy results per length
     * @return         The optimal length (maximizing entropy)
     */
    static int selectBestLength(const std::vector<EntropyResult>& results);

    /**
     * Select TOP-K lengths for multi-length aggregation.
     * Combines MAWs from multiple high-entropy lengths to increase
     * the number of informative characters — critical for:
     *   - Short sequences (<1kb) where single-length MAWs are sparse
     *   - Highly similar species where single-length gives few variable chars
     *
     * @param results  Entropy results per length
     * @param topK     Number of top lengths to select (default: 3)
     * @param minChars Minimum variable characters required (default: 10)
     * @return         Vector of selected lengths, sorted
     */
    static std::vector<int> selectTopLengths(const std::vector<EntropyResult>& results,
                                              int topK = 3, int minChars = 10);

    /**
     * Compute adaptive length range based on sequence properties.
     * For short sequences: uses smaller lengths (2-6)
     * For long sequences: uses wider range (3-12)
     * For highly similar species: expands range to capture divergence
     *
     * @param avgSeqLen  Average sequence length
     * @param numSpecies Number of species
     * @param lmin       Output: recommended minimum length
     * @param lmax       Output: recommended maximum length
     */
    static void computeAdaptiveRange(int avgSeqLen, int numSpecies,
                                      int& lmin, int& lmax);
};

} // namespace mlmaws

#endif // ENTROPY_SELECTOR_H
