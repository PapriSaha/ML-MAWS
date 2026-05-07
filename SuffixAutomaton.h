/**
 * ML-MAWS: Maximum Likelihood Phylogeny Estimation Using
 *          Minimal Absent Word Characters
 * 
 * SuffixAutomaton.h - Suffix Automaton for MAW extraction
 * 
 * Based on the suffix automaton approach from Ehsan et al. (WALCOM 2025)
 * "An Efficient Implementation of CD-MAWS using Suffix Automata"
 * 
 * Authors: [Your Name]
 * License: MIT
 */

#ifndef SUFFIX_AUTOMATON_H
#define SUFFIX_AUTOMATON_H

#include <vector>
#include <string>
#include <cstring>

namespace mlmaws {

// DNA alphabet: A=0, C=1, G=2, T=3
const int SIGMA = 4;
const std::string ALPHABET = "ACGT";

/**
 * Vertex (state) in the Suffix Automaton.
 */
struct SAVertex {
    int link;       // suffix link
    int len;        // longest string length reaching this state
    int next[SIGMA]; 

    SAVertex() : link(-1), len(0) {
        memset(next, -1, sizeof(next));
    }
};

/**
 * Suffix Automaton data structure.
 * Supports online construction in O(n) time for a string of length n.
 */
class SuffixAutomaton {
public:
    SuffixAutomaton();
    ~SuffixAutomaton() = default;

    /** Reset the automaton for a new string. */
    void init();

    /** Add a single character (A/C/G/T) to the automaton. */
    void addChar(char c);

    /** Build the automaton for an entire DNA string. */
    void build(const std::string& s);

    /** Access states directly. */
    const SAVertex& state(int idx) const { return sa_[idx]; }
    int size() const { return sz_; }

    /** Map DNA character to integer 0..3. */
    static int charToInt(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }

    /** Map integer 0..3 to DNA character. */
    static char intToChar(int i) {
        return ALPHABET[i];
    }

private:
    std::vector<SAVertex> sa_;
    int last_;
    int sz_;

    void ensureCapacity(int idx);
};

} // namespace mlmaws

#endif // SUFFIX_AUTOMATON_H
