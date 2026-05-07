/**
 * MAWExtractor.h - Minimal Absent Word extraction using Suffix Automaton
 * 
 * Extracts MAWs within a specified length range [minK, maxK].
 * Supports strand-aware filtering (intersection of forward and reverse complement MAWs).
 */

#ifndef MAW_EXTRACTOR_H
#define MAW_EXTRACTOR_H

#include "SuffixAutomaton.h"
#include <vector>
#include <string>
#include <set>

namespace mlmaws {

/**
 * Extracts Minimal Absent Words from a DNA sequence using a Suffix Automaton.
 */
class MAWExtractor {
public:
    MAWExtractor(int minK = 2, int maxK = 12);
    ~MAWExtractor() = default;

    /** Set length range for MAW extraction. */
    void setLengthRange(int minK, int maxK);

    /**
     * Extract MAWs from a DNA sequence.
     * Returns sorted vector of MAW strings.
     */
    std::vector<std::string> extract(const std::string& sequence);

    /**
     * Extract strand-aware MAWs: MAW must be absent from
     * BOTH the forward strand and its reverse complement.
     * Returns sorted vector of MAW strings.
     */
    std::vector<std::string> extractStrandAware(const std::string& sequence);

    /**
     * Extract MAWs restricted to a specific length.
     * Returns sorted vector of MAW strings of exactly that length.
     */
    std::vector<std::string> extractByLength(const std::string& sequence, int length);

    /** Get reverse complement of a DNA string. */
    static std::string reverseComplement(const std::string& seq);

private:
    int minK_;
    int maxK_;
    SuffixAutomaton sa_;

    // DFS state
    char buffer_[64]; // max MAW length buffer
    int bufLen_;
    std::vector<std::string>* currentResult_;

    /** DFS traversal of suffix automaton to enumerate MAWs. */
    void dfs(int node);
};

} // namespace mlmaws

#endif // MAW_EXTRACTOR_H
