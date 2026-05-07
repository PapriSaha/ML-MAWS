/**
 * MAWExtractor.cpp - Minimal Absent Word extraction implementation
 */

#include "MAWExtractor.h"
#include <algorithm>
#include <iostream>

namespace mlmaws {

MAWExtractor::MAWExtractor(int minK, int maxK)
    : minK_(minK), maxK_(maxK), bufLen_(0), currentResult_(nullptr) {}

void MAWExtractor::setLengthRange(int minK, int maxK) {
    minK_ = minK;
    maxK_ = maxK;
}

std::string MAWExtractor::reverseComplement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    for (int i = (int)seq.size() - 1, j = 0; i >= 0; --i, ++j) {
        switch (seq[i]) {
            case 'A': case 'a': rc[j] = 'T'; break;
            case 'C': case 'c': rc[j] = 'G'; break;
            case 'G': case 'g': rc[j] = 'C'; break;
            case 'T': case 't': rc[j] = 'A'; break;
            default: rc[j] = 'N'; break;
        }
    }
    return rc;
}

void MAWExtractor::dfs(int node) {
    if (bufLen_ >= maxK_) return;
    if (bufLen_ >= 62) return; // buffer overflow guard

    const SAVertex& st = sa_.state(node);
    int link = st.link;

    for (int i = 0; i < SIGMA; i++) {
        int u = st.next[i];
        if (u == -1) {
            // Character i leads to no state => potential MAW
            if (bufLen_ + 1 < minK_) continue;

            // Check: the suffix link's transition must exist and
            // the suffix link's length must be >= current depth
            if (link >= 0) {
                const SAVertex& linkSt = sa_.state(link);
                int v = linkSt.next[i];
                if (v != -1 && linkSt.len + 1 >= bufLen_) {
                    buffer_[bufLen_] = SuffixAutomaton::intToChar(i);
                    buffer_[bufLen_ + 1] = '\0';
                    currentResult_->emplace_back(buffer_, bufLen_ + 1);
                }
            }
        } else {
            // Character i has a transition => go deeper
            buffer_[bufLen_++] = SuffixAutomaton::intToChar(i);
            dfs(u);
            bufLen_--;
        }
    }
}

std::vector<std::string> MAWExtractor::extract(const std::string& sequence) {
    std::vector<std::string> result;

    sa_.build(sequence);
    bufLen_ = 0;
    currentResult_ = &result;
    dfs(0);
    currentResult_ = nullptr;

    // DFS does NOT produce sorted output (SA DAG is not lexicographic)
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<std::string> MAWExtractor::extractStrandAware(const std::string& sequence) {
    // Get MAWs from forward strand
    std::vector<std::string> fwdMAWs = extract(sequence);

    // Get MAWs from reverse complement
    std::string rc = reverseComplement(sequence);
    std::vector<std::string> rcMAWs = extract(rc);

    // Intersection: keep only MAWs present in BOTH sets
    // Both vectors are already sorted
    std::vector<std::string> result;
    std::set_intersection(
        fwdMAWs.begin(), fwdMAWs.end(),
        rcMAWs.begin(), rcMAWs.end(),
        std::back_inserter(result)
    );

    return result;
}

std::vector<std::string> MAWExtractor::extractByLength(const std::string& sequence, int length) {
    int origMin = minK_, origMax = maxK_;
    minK_ = length;
    maxK_ = length;

    std::vector<std::string> all = extract(sequence);

    minK_ = origMin;
    maxK_ = origMax;

    // Filter to exact length (should already be exact, but just in case)
    std::vector<std::string> result;
    for (auto& w : all) {
        if ((int)w.size() == length) {
            result.push_back(std::move(w));
        }
    }
    return result;
}

} // namespace mlmaws
