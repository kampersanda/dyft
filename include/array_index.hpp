#pragma once

#include "dyft_factors.hpp"
#include "dyft_interface.hpp"
#include "sparse_table.hpp"
#include "vcode_array.hpp"

namespace dyft {

template <int N>
class array_index : public dyft_interface<N> {
  public:
    using vint_type = typename vcode_tools<N>::vint_type;

    static constexpr int LEN = N;
    static constexpr int NOT_FOUND = std::numeric_limits<int>::max();

    static int get_chunks_size(int) {
        return LEN;
    }

  private:
    const vcode_array<N>* m_database = nullptr;

    const int m_sigma = 0;
    const int m_radius = 0;
    const int m_splitthr = 0;
    const int* m_splitthrs_ptr = nullptr;

    // For computing trie cost
    const double* m_infacts_ptr = nullptr;
    const double* m_outfacts_ptr = nullptr;

    // Trie
    std::vector<int> m_fanouts;
    sparse_table m_idseqs;
    uint32_t m_ids = 0;  // num ids

    const int m_depth_beg = 0;
    const int m_depth_end = N;

    // double m_bf_cost = 0.0;
    double m_trie_cost = 0.0;
    double m_in_weight = 1.0;

    // For query processing
    struct state_type {
        int npos;
        int depth;
        int dist;
    };
    std::vector<state_type> m_states;

    // statistics
    size_t m_split_count = 0;

  public:
    explicit array_index(const vcode_array<N>* database, int radius, int splitthr, double in_weight)
        : array_index(database, radius, 0, N, splitthr, in_weight) {}

    explicit array_index(const vcode_array<N>* database, int radius, int depth_beg, int depth_end,  //
                         int splitthr, double in_weight)
        : m_database(database), m_sigma(1 << m_database->get_bits()), m_radius(radius), m_splitthr(splitthr),
          m_splitthrs_ptr(m_splitthr <= 0 ? dyft_factors::get_splitthrs(m_database->get_bits(), m_radius) : nullptr),
          m_infacts_ptr(m_splitthr <= 0 ? dyft_factors::get_infacts(m_database->get_bits(), m_radius) : nullptr),
          m_outfacts_ptr(m_splitthr <= 0 ? dyft_factors::get_outfacts(m_database->get_bits(), m_radius) : nullptr),
          m_fanouts(m_sigma, NOT_FOUND), m_depth_beg(depth_beg), m_depth_end(depth_end), m_in_weight(in_weight) {
        m_states.reserve(1 << 10);
    }

    uint32_t append() override {
        ABORT_IF_LE(m_database->get_size(), m_ids);

        const uint32_t new_id = m_ids++;
        const vint_type* vcode = m_database->access(new_id);

        int npos = 0, depth = m_depth_beg;

        while (true) {
            const int c = vcode_tools<N>::get_int(vcode, depth++, get_bits());
            const int cpos = npos + c;

            if (m_fanouts[cpos] == NOT_FOUND) {
                npos = cpos;
                m_fanouts[npos] = -1 * int(m_idseqs.size());
                m_idseqs.push_back();
                break;
            }

            if (m_fanouts[cpos] <= 0) {  // Leaf?
                npos = cpos;
                break;
            }

            npos = m_fanouts[cpos] * m_sigma;
        }

        const int lpos = -1 * m_fanouts[npos];
        m_idseqs.insert(lpos, new_id);

        if (m_outfacts_ptr) {
            m_trie_cost += m_outfacts_ptr[depth - m_depth_beg - 1];
        }

        if (depth == m_depth_end) {
            return new_id;
        }

        const int cnt = int(m_idseqs.group_size(lpos));

        if (m_splitthrs_ptr == nullptr) {
            if (m_splitthr <= cnt) {
                split_node(npos, depth);
            }
        } else {
            if (m_splitthrs_ptr[depth - 1] * m_in_weight <= cnt) {
                split_node(npos, depth);
            }
        }
        return new_id;
    }

    void range_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) override {
        return selects_ls() ? ls_search(vcode, fn) : trie_search(vcode, fn);
    }

    void trie_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) {
        m_states.clear();
        m_states.push_back(state_type{0, m_depth_beg, 0});

        int code[N];
        for (int i = m_depth_beg; i < m_depth_end; i++) {
            code[i] = vcode_tools<N>::get_int(vcode, i, get_bits());
        }

        while (!m_states.empty()) {
            state_type s = m_states.back();
            m_states.pop_back();

            if (s.dist < m_radius) {
                const int c = code[s.depth];

                for (int i = 0; i < m_sigma; ++i) {
                    const int cpos = s.npos + i;

                    if (m_fanouts[cpos] == NOT_FOUND) {
                        continue;
                    }

                    if (m_fanouts[cpos] <= 0) {  // leaf
                        const int lpos = -1 * m_fanouts[cpos];
                        auto [bptr, eptr] = m_idseqs.access(lpos);
                        for (auto it = bptr; it != eptr; ++it) {
                            fn(*it);
                        }
                    } else {  // internal
                        if (i == c) {
                            m_states.push_back(state_type{m_fanouts[cpos] * m_sigma, s.depth + 1, s.dist});
                        } else {
                            m_states.push_back(state_type{m_fanouts[cpos] * m_sigma, s.depth + 1, s.dist + 1});
                        }
                    }
                }
            } else {  // exact match
                while (true) {
                    const int c = code[s.depth];
                    const int cpos = s.npos + c;

                    if (m_fanouts[cpos] == NOT_FOUND) {
                        break;
                    }

                    if (m_fanouts[cpos] <= 0) {  // leaf
                        const int lpos = -1 * m_fanouts[cpos];
                        auto [bptr, eptr] = m_idseqs.access(lpos);
                        for (auto it = bptr; it != eptr; ++it) {
                            fn(*it);
                        }
                        break;
                    }

                    // internal
                    s.npos = m_fanouts[cpos] * m_sigma;
                    s.depth += 1;
                }
            }
        }
    }

    void ls_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) {
        for (uint32_t id = 0; id < get_size(); id++) {
            fn(id);
        }
    }

    uint32_t get_size() const override {
        return m_ids;
    }
    uint32_t get_leaves() const override {
        return m_idseqs.size();
    }
    int get_bits() const override {
        return m_database->get_bits();
    }

    double get_ls_cost() const override {
        return static_cast<double>(m_ids) * get_bits();
    }
    double get_trie_cost() const override {
        return m_trie_cost;
    }

    bool selects_ls() const override {
        return m_infacts_ptr and get_ls_cost() < get_trie_cost();
    }

    size_t get_split_count() override {
        return m_split_count;
    }

    void innode_stats() override {
        tfm::printfln("innode_stats: %d", m_fanouts.size() / m_sigma);
    }
    void population_stats() override {}

  private:
    void split_node(const int npos, const int depth) {
        ABORT_IF_LT(0, m_fanouts[npos]);
        ABORT_IF_LT(depth, m_depth_beg + 1);

        m_split_count += 1;

        if (m_fanouts.size() == m_fanouts.capacity()) {
            m_fanouts.reserve(m_fanouts.size() * 2);
        }

        const int new_npos = static_cast<int>(m_fanouts.size());
        m_fanouts.resize(new_npos + m_sigma, NOT_FOUND);

        int lpos = -1 * m_fanouts[npos];
        ABORT_IF_LE(m_idseqs.size(), size_t(lpos));

        // update cost
        if (m_infacts_ptr) {
            const int _depth = depth - m_depth_beg - 1;
            const int g_size = m_idseqs.group_size(lpos);

            m_trie_cost -= g_size * m_outfacts_ptr[_depth];
            m_trie_cost += (m_infacts_ptr[_depth] * m_in_weight) + (g_size * m_outfacts_ptr[_depth + 1]);
        }

        std::vector<std::vector<uint32_t> > idbufs(m_sigma);
        {
            const auto ids = m_idseqs.extract(lpos);
            for (uint32_t id : ids) {
                const vint_type* vcode = m_database->access(id);
                const int c = vcode_tools<N>::get_int(vcode, depth, get_bits());
                idbufs[c].push_back(id);
            }
        }

        for (int c = 0; c < m_sigma; c++) {
            if (idbufs[c].empty()) {
                continue;
            }

            const int new_cpos = new_npos + c;

            if (lpos != NOT_FOUND) {
                // the first new idseq
                m_fanouts[new_cpos] = -1 * lpos;
                m_idseqs.insert(lpos, idbufs[c]);
                lpos = NOT_FOUND;
            } else {
                m_fanouts[new_cpos] = -1 * int(m_idseqs.size());
                m_idseqs.push_back(idbufs[c]);
            }
        }

        m_fanouts[npos] = new_npos / m_sigma;
    }
};

}  // namespace dyft
