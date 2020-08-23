#pragma once

#include "dyft_factors.hpp"
#include "dyft_interface.hpp"
#include "mart_array_dense.hpp"
#include "mart_array_full.hpp"
#include "mart_array_sparse.hpp"
#include "sparse_table.hpp"
#include "vcode_array.hpp"

namespace dyft {

template <int N>
class art_index : public dyft_interface<N> {
  public:
    using vint_type = typename vcode_tools<N>::vint_type;

    using array_4_type = mart_array_sparse<4, MART_NODE_4>;
    using array_16_type = mart_array_sparse<16, MART_NODE_16>;
    using array_48_type = mart_array_dense<48, MART_NODE_64>;  // instead of MART_NODE_48
    using array_256_type = mart_array_full<MART_NODE_256>;

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

    // ART
    array_4_type m_array_4;
    array_16_type m_array_16;
    array_48_type m_array_48;
    array_256_type m_array_256;

    mart_pointer m_rootptr;

    // std::vector<int> m_fanouts;
    sparse_table m_idseqs;
    uint32_t m_ids = 0;  // num ids

    const int m_depth_beg = 0;
    const int m_depth_end = N;

    // double m_bf_cost = 0.0;
    double m_trie_cost = 0.0;
    double m_in_weight = 1.0;

    // For query processing
    struct state_type {
        mart_pointer nptr;
        int depth;
        int dist;
    };
    std::vector<state_type> m_states;

    array_4_type::searcher m_s4;
    array_16_type::searcher m_s16;
    array_48_type::searcher m_s48;
    array_256_type::searcher m_s256;
    std::vector<mart_edge> m_edges;

    // statistics
    size_t m_split_count = 0;

  public:
    explicit art_index(const vcode_array<N>* database, int radius, int splitthr, double in_weight)
        : art_index(database, radius, 0, N, splitthr, in_weight) {}

    explicit art_index(const vcode_array<N>* database, int radius, int depth_beg, int depth_end,  //
                       int splitthr, double in_weight)
        : m_database(database), m_sigma(1 << m_database->get_bits()), m_radius(radius), m_splitthr(splitthr),
          m_splitthrs_ptr(m_splitthr <= 0 ? dyft_factors::get_splitthrs(m_database->get_bits(), m_radius) : nullptr),
          m_infacts_ptr(m_splitthr <= 0 ? dyft_factors::get_infacts(m_database->get_bits(), m_radius) : nullptr),
          m_outfacts_ptr(m_splitthr <= 0 ? dyft_factors::get_outfacts(m_database->get_bits(), m_radius) : nullptr),
          //   m_fanouts(m_sigma, NOT_FOUND),
          m_depth_beg(depth_beg), m_depth_end(depth_end), m_in_weight(in_weight) {
        m_rootptr = m_array_256.make_node();  // make root
        m_states.reserve(1 << 10);
        m_edges.reserve(256);
    }

    uint32_t append() override {
        ABORT_IF_LE(m_database->get_size(), m_ids);

        const uint32_t new_id = m_ids++;
        const vint_type* vcode = m_database->access(new_id);

        // int npos = 0, depth = m_depth_beg;
        int depth = m_depth_beg;
        mart_cursor mc{UINT32_MAX, make_mart_nullptr(), m_rootptr};

        const uint32_t new_lpos = m_idseqs.size();
        const mart_pointer new_lptr{new_lpos, MART_LEAF};

        while (true) {
            DEBUG_ABORT_IF_LE(m_depth_end, depth);
            const int c = vcode_tools<N>::get_int(vcode, depth++, get_bits());

            mart_insert_flags iflag;

            switch (mc.nptr.ntype) {
                case MART_NODE_4: {
                    iflag = m_array_4.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_16: {
                    iflag = m_array_16.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_64: {  // instead of MART_NODE_48
                    iflag = m_array_48.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_256: {
                    iflag = m_array_256.insert_ptr(mc, c, new_lptr);
                    break;
                }
                default: {
                    ABORT_IF(true);  // should not come
                    break;
                }
            }

            if (iflag == MART_INSERTED) {
                ABORT_IF_NE(mc.nptr.nid, new_lpos);
                ABORT_IF_NE(mc.nptr.ntype, MART_LEAF);
                m_idseqs.push_back();
                break;
            }

            if (mc.nptr.ntype == MART_LEAF) {
                break;
            }
        }

        const uint32_t lpos = mc.nptr.nid;
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
                split_node(mc, depth);
            }
        } else {
            if (m_splitthrs_ptr[depth - 1] * m_in_weight <= cnt) {
                split_node(mc, depth);
            }
        }
        return new_id;
    }

    void range_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) override {
        return selects_ls() ? ls_search(vcode, fn) : trie_search(vcode, fn);
    }

    void trie_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) {
        m_states.clear();
        m_states.push_back(state_type{m_rootptr, m_depth_beg, 0});

        int code[N];
        for (int i = m_depth_beg; i < m_depth_end; i++) {
            code[i] = vcode_tools<N>::get_int(vcode, i, get_bits());
        }

        while (!m_states.empty()) {
            state_type s = m_states.back();
            m_states.pop_back();

            if (s.dist < m_radius) {
                const uint8_t c = code[s.depth];
                m_edges.clear();

                switch (s.nptr.ntype) {
                    case MART_NODE_4: {
                        m_s4.set(&m_array_4, s.nptr);
                        m_s4.scan(m_edges);
                        break;
                    }
                    case MART_NODE_16: {
                        m_s16.set(&m_array_16, s.nptr);
                        m_s16.scan(m_edges);
                        break;
                    }
                    case MART_NODE_64: {
                        m_s48.set(&m_array_48, s.nptr);
                        m_s48.scan(m_edges);
                        break;
                    }
                    case MART_NODE_256: {
                        m_s256.set(&m_array_256, s.nptr);
                        m_s256.scan(m_edges);
                        break;
                    }
                    default: {
                        ABORT_IF(true);  // should not come
                        break;
                    }
                }

                for (const mart_edge& e : m_edges) {
                    if (e.nptr.ntype == MART_LEAF) {  // leaf
                        const uint32_t lpos = e.nptr.nid;
                        auto [bptr, eptr] = m_idseqs.access(lpos);
                        for (auto it = bptr; it != eptr; ++it) {
                            fn(*it);
                        }
                    } else {  // internal
                        if (e.c == c) {
                            m_states.push_back(state_type{e.nptr, s.depth + 1, s.dist});
                        } else {
                            m_states.push_back(state_type{e.nptr, s.depth + 1, s.dist + 1});
                        }
                    }
                }
            } else {  // exact match
                while (true) {
                    const uint8_t c = code[s.depth];

                    switch (s.nptr.ntype) {
                        case MART_NODE_4: {
                            s.nptr = m_array_4.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_16: {
                            s.nptr = m_array_16.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_64: {
                            s.nptr = m_array_48.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_256: {
                            s.nptr = m_array_256.find_child(s.nptr, c);
                            break;
                        }
                        default: {
                            ABORT_IF(true);  // should not come
                            break;
                        }
                    }

                    if (s.nptr.nid == MART_NILID) {
                        break;
                    }

                    if (s.nptr.ntype == MART_LEAF) {
                        const uint32_t lpos = s.nptr.nid;
                        auto [bptr, eptr] = m_idseqs.access(lpos);
                        for (auto it = bptr; it != eptr; ++it) {
                            fn(*it);
                        }
                        break;
                    }

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
        tfm::printfln("## innode_stats ##");
        tfm::printfln("- num node 4:   %d (%d)", m_array_4.num_nodes(), m_array_4.num_emp_nodes());
        tfm::printfln("- num node 16:  %d (%d)", m_array_16.num_nodes(), m_array_16.num_emp_nodes());
        tfm::printfln("- num node 48:  %d (%d)", m_array_48.num_nodes(), m_array_48.num_emp_nodes());
        tfm::printfln("- num node 256: %d (%d)", m_array_256.num_nodes(), m_array_256.num_emp_nodes());
    }
    void population_stats() override {}

  private:
    mart_insert_flags expand(mart_cursor& mc, uint8_t c, mart_pointer new_ptr) {
        m_edges.clear();

        switch (mc.nptr.ntype) {
            case MART_NODE_4: {
                m_array_4.extract(mc, m_edges);
                mc.nptr = m_array_16.make_node(m_edges);
                update_srcptr(mc);
                return m_array_16.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_16: {
                m_array_16.extract(mc, m_edges);
                mc.nptr = m_array_48.make_node(m_edges);
                update_srcptr(mc);
                return m_array_48.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_64: {  // instead of MART_NODE_48
                m_array_48.extract(mc, m_edges);
                mc.nptr = m_array_256.make_node(m_edges);
                update_srcptr(mc);
                return m_array_256.insert_ptr(mc, c, new_ptr);
            }
            default: {
                ABORT_IF(true);  // should not come
                break;
            }
        }
    }

    void update_srcptr(const mart_cursor& mc) {
        switch (mc.pptr.ntype) {
            case MART_NODE_4: {
                m_array_4.update_srcptr(mc);
                break;
            }
            case MART_NODE_16: {
                m_array_16.update_srcptr(mc);
                break;
            }
            case MART_NODE_64: {  // instead of MART_NODE_48
                m_array_48.update_srcptr(mc);
                break;
            }
            case MART_NODE_256: {
                m_array_256.update_srcptr(mc);
                break;
            }
            default: {
                ABORT_IF(true);  // should not come
                break;
            }
        }
    }

    void split_node(mart_cursor& mc, const int depth) {
        ABORT_IF_LE(depth, m_depth_beg);
        ABORT_IF_NE(mc.nptr.ntype, MART_LEAF);

        m_split_count += 1;

        const uint32_t lpos = mc.nptr.nid;
        ABORT_IF_LE(m_idseqs.size(), lpos);

        // update cost
        if (m_infacts_ptr) {
            const int _depth = depth - m_depth_beg - 1;
            const int g_size = m_idseqs.group_size(lpos);

            m_trie_cost -= g_size * m_outfacts_ptr[_depth];
            m_trie_cost += (m_infacts_ptr[_depth] * m_in_weight) + (g_size * m_outfacts_ptr[_depth + 1]);
        }

        std::vector<std::vector<uint32_t>> idbufs(m_sigma);
        {
            const auto ids = m_idseqs.extract(lpos);
            for (uint32_t id : ids) {
                const vint_type* vcode = m_database->access(id);
                const int c = vcode_tools<N>::get_int(vcode, depth, get_bits());
                idbufs[c].push_back(id);
            }
        }

        m_edges.clear();

        for (int c = 0; c < m_sigma; c++) {
            if (idbufs[c].empty()) {
                continue;
            }
            if (m_edges.empty()) {
                const mart_pointer new_lptr = mart_pointer{lpos, MART_LEAF};
                m_edges.push_back(mart_edge{uint8_t(c), new_lptr});
                m_idseqs.insert(lpos, idbufs[c]);
            } else {
                const mart_pointer new_lptr = mart_pointer{uint32_t(m_idseqs.size()), MART_LEAF};
                m_edges.push_back(mart_edge{uint8_t(c), new_lptr});
                m_idseqs.push_back(idbufs[c]);
            }
        }

        if (m_edges.size() <= 4) {
            mc.nptr = m_array_4.make_node(m_edges);
        } else if (m_edges.size() <= 16) {
            mc.nptr = m_array_16.make_node(m_edges);
        } else if (m_edges.size() <= 48) {
            mc.nptr = m_array_48.make_node(m_edges);
        } else {
            mc.nptr = m_array_256.make_node(m_edges);
        }

        update_srcptr(mc);
    }
};

}  // namespace dyft
