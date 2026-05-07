/**
 * SuffixAutomaton.cpp - Suffix Automaton implementation
 */

#include "SuffixAutomaton.h"

namespace mlmaws {

SuffixAutomaton::SuffixAutomaton() : last_(0), sz_(1) {
    sa_.resize(16);
}

void SuffixAutomaton::ensureCapacity(int idx) {
    while (idx >= (int)sa_.size()) {
        sa_.resize(sa_.size() * 2);
    }
}

void SuffixAutomaton::init() {
    for (int i = 0; i < sz_; i++) {
        sa_[i] = SAVertex();
    }
    last_ = 0;
    sz_ = 1;
}

void SuffixAutomaton::addChar(char c) {
    int ci = charToInt(c);
    if (ci < 0) return; // skip non-DNA characters

    int cur = sz_++;
    ensureCapacity(cur);
    sa_[cur].len = sa_[last_].len + 1;

    int u = last_;
    while (u != -1 && sa_[u].next[ci] == -1) {
        sa_[u].next[ci] = cur;
        u = sa_[u].link;
    }

    if (u == -1) {
        sa_[cur].link = 0;
    } else {
        int v = sa_[u].next[ci];
        if (sa_[u].len + 1 == sa_[v].len) {
            sa_[cur].link = v;
        } else {
            int nw = sz_++;
            ensureCapacity(nw);
            sa_[nw].link = sa_[v].link;
            sa_[nw].len = sa_[u].len + 1;
            memcpy(sa_[nw].next, sa_[v].next, sizeof(sa_[v].next));
            while (u != -1 && sa_[u].next[ci] == v) {
                sa_[u].next[ci] = nw;
                u = sa_[u].link;
            }
            sa_[cur].link = sa_[v].link = nw;
        }
    }
    last_ = cur;
}

void SuffixAutomaton::build(const std::string& s) {
    init();
    ensureCapacity((int)s.size() * 2 + 5);
    for (char c : s) {
        addChar(c);
    }
}

} // namespace mlmaws
