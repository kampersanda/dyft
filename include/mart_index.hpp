#pragma once

#include "dyft_factors.hpp"
#include "dyft_interface.hpp"
#include "ham_tables.hpp"
#include "mart_array_dense.hpp"
#include "mart_array_full.hpp"
#include "mart_array_sparse.hpp"
#include "misc.hpp"
#include "sparse_table.hpp"
#include "vcode_array.hpp"

namespace dyft {

template <int N>
class mart_index : public dyft_interface<N> {
  public:
    using vint_type = typename vcode_tools<N>::vint_type;

#ifdef MART_SPACE_EFFICIENT
    using array_2_type = mart_array_sparse<2, MART_NODE_2>;
    using array_4_type = mart_array_sparse<4, MART_NODE_4>;
    using array_8_type = mart_array_sparse<8, MART_NODE_8>;
    using array_16_type = mart_array_sparse<16, MART_NODE_16>;
    using array_32_type = mart_array_sparse<32, MART_NODE_32>;
    using array_64_type = mart_array_dense<64, MART_NODE_64>;
    using array_128_type = mart_array_dense<128, MART_NODE_128>;
    using array_256_type = mart_array_full<MART_NODE_256>;
#else
    using array_4_type = mart_array_sparse<4, MART_NODE_4>;
    using array_32_type = mart_array_sparse<32, MART_NODE_32>;
    using array_256_type = mart_array_full<MART_NODE_256>;
#endif

    static constexpr int LEN = N;
    static constexpr int NOT_FOUND = std::numeric_limits<int>::max();

    static int get_chunks_size(int bits) {
        return ceil_div(N, 8 / bits);
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

    // MART
#ifdef MART_SPACE_EFFICIENT
    array_2_type m_array_2;
    array_4_type m_array_4;
    array_8_type m_array_8;
    array_16_type m_array_16;
    array_32_type m_array_32;
    array_64_type m_array_64;
    array_128_type m_array_128;
    array_256_type m_array_256;
#else
    array_4_type m_array_4;
    array_32_type m_array_32;
    array_256_type m_array_256;
#endif

    mart_pointer m_rootptr;

    sparse_table m_idseqs;  // positing lists
    uint32_t m_ids = 0;  // num ids

    // For packing
    uint8_t m_chunks[N];  // of length CHUNKS
    const int m_ints_per_chunk;  // num ints per chunk
    const int m_chunks_size;

    const int m_bpos_beg;
    const int m_bpos_end;

    double m_trie_cost;
    double m_in_weight;

    // For query processing
    struct state_type {
        mart_pointer nptr;
        int bpos;
        int dist;
    };
    std::vector<state_type> m_states;

#ifdef MART_SPACE_EFFICIENT
    array_2_type::searcher m_s2;
    array_4_type::searcher m_s4;
    array_8_type::searcher m_s8;
    array_16_type::searcher m_s16;
    array_32_type::searcher m_s32;
    array_64_type::searcher m_s64;
    array_128_type::searcher m_s128;
    array_256_type::searcher m_s256;
#else
    array_4_type::searcher m_s4;
    array_32_type::searcher m_s32;
    array_256_type::searcher m_s256;
#endif
    std::vector<mart_edge> m_edges;

    // statistics
    size_t m_split_count = 0;

  public:
    explicit mart_index(const vcode_array<N>* database, int radius, int splitthr, double in_weight)
        : mart_index(database, radius, 0, get_chunks_size(database->get_bits()), splitthr, in_weight) {}

    explicit mart_index(const vcode_array<N>* database, int radius, int bpos_beg, int bpos_end,  //
                        int splitthr, double in_weight)
        : m_database(database), m_sigma(1 << m_database->get_bits()), m_radius(radius), m_splitthr(splitthr),
          m_splitthrs_ptr(m_splitthr <= 0 ? dyft_factors::get_splitthrs(m_database->get_bits(), m_radius) : nullptr),
          m_infacts_ptr(m_splitthr <= 0 ? dyft_factors::get_infacts(m_database->get_bits(), m_radius) : nullptr),
          m_outfacts_ptr(m_splitthr <= 0 ? dyft_factors::get_outfacts(m_database->get_bits(), m_radius) : nullptr),
          m_ints_per_chunk(8 / m_database->get_bits()), m_chunks_size(ceil_div(N, m_ints_per_chunk)),
          m_bpos_beg(bpos_beg), m_bpos_end(bpos_end), m_trie_cost(0.0), m_in_weight(in_weight) {
        m_rootptr = m_array_256.make_node();  // make root
        m_states.reserve(1 << 10);
        m_edges.reserve(256);
    }

    uint32_t append() override {
        ABORT_IF_LE(m_database->get_size(), m_ids);

        const uint32_t new_id = m_ids++;
        const vint_type* vcode = m_database->access(new_id);

        to_chunks(vcode);

        int bpos = m_bpos_beg;  // chunk position
        mart_cursor mc{UINT32_MAX, make_mart_nullptr(), m_rootptr};

        const uint32_t new_lpos = m_idseqs.size();
        const mart_pointer new_lptr{new_lpos, MART_LEAF};

        while (true) {
            DEBUG_ABORT_IF_LE(m_bpos_end, bpos);
            const uint8_t c = m_chunks[bpos++];

            mart_insert_flags iflag;

#ifdef MART_SPACE_EFFICIENT
            switch (mc.nptr.ntype) {
                case MART_NODE_2: {
                    iflag = m_array_2.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_4: {
                    iflag = m_array_4.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_8: {
                    iflag = m_array_8.insert_ptr(mc, c, new_lptr);
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
                case MART_NODE_32: {
                    iflag = m_array_32.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_64: {
                    iflag = m_array_64.insert_ptr(mc, c, new_lptr);
                    if (iflag == MART_NEEDED_TO_EXPAND) {
                        iflag = expand(mc, c, new_lptr);
                        ABORT_IF_NE(iflag, MART_INSERTED);
                    }
                    break;
                }
                case MART_NODE_128: {
                    iflag = m_array_128.insert_ptr(mc, c, new_lptr);
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
#else
            if (mc.nptr.ntype == MART_NODE_4) {
                iflag = m_array_4.insert_ptr(mc, c, new_lptr);
                if (iflag == MART_NEEDED_TO_EXPAND) {
                    iflag = expand(mc, c, new_lptr);
                    ABORT_IF_NE(iflag, MART_INSERTED);
                }
            } else if (mc.nptr.ntype == MART_NODE_32) {
                iflag = m_array_32.insert_ptr(mc, c, new_lptr);
                if (iflag == MART_NEEDED_TO_EXPAND) {
                    iflag = expand(mc, c, new_lptr);
                    ABORT_IF_NE(iflag, MART_INSERTED);
                }
            } else if (mc.nptr.ntype == MART_NODE_256) {
                iflag = m_array_256.insert_ptr(mc, c, new_lptr);
            } else {
                ABORT_IF(true);  // should not come
            }
#endif

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
            const int _depth = (bpos - m_bpos_beg - 1) * m_ints_per_chunk;
            m_trie_cost += m_outfacts_ptr[_depth];
            // m_trie_cost += m_outfacts_ptr[bpos - m_bpos_beg - 1];
        }

        if (bpos == m_bpos_end) {
            return new_id;
        }

        const int cnt = int(m_idseqs.group_size(lpos));

        if (m_splitthrs_ptr == nullptr) {
            if (m_splitthr <= cnt) {
                split_node(mc, bpos);
            }
        } else {
            const int _depth = (bpos - m_bpos_beg - 1) * m_ints_per_chunk;
            if (m_splitthrs_ptr[_depth] * m_in_weight <= cnt) {
                split_node(mc, bpos);
            }
            // if (m_splitthrs_ptr[_depth] <= cnt) {
            //     split_node(mc, bpos);
            // }
        }
        return new_id;
    }

    void range_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) override {
        return selects_ls() ? ls_search(vcode, fn) : trie_search(vcode, fn);
    }

    void trie_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) {
        m_states.clear();
        m_states.push_back(state_type{m_rootptr, m_bpos_beg, 0});

        to_chunks(vcode);
        const int* bpt = ham_tables::get_bpt(get_bits());

        while (!m_states.empty()) {
            state_type s = m_states.back();
            m_states.pop_back();

            DEBUG_ABORT_IF_LT(m_radius, s.dist);
            DEBUG_ABORT_IF_LE(m_bpos_end, s.bpos);

            if (s.dist < m_radius) {
                const uint8_t c = m_chunks[s.bpos];
                const uint8_t* lut = ham_tables::get_lut(get_bits(), c);
                const uint8_t* hdt = ham_tables::get_hdt(get_bits(), c);

                m_edges.clear();
                const int r = std::min(m_radius - s.dist, m_ints_per_chunk);

#ifdef MART_SPACE_EFFICIENT
                switch (s.nptr.ntype) {
                    case MART_NODE_2: {
                        m_s2.set(&m_array_2, s.nptr);
                        m_s2.scan(m_edges);
                        break;
                    }
                    case MART_NODE_4: {
                        m_s4.set(&m_array_4, s.nptr);
                        m_s4.scan(m_edges);
                        break;
                    }
                    case MART_NODE_8: {
                        m_s8.set(&m_array_8, s.nptr);
                        m_s8.scan(m_edges);
                        break;
                    }
                    case MART_NODE_16: {
                        m_s16.set(&m_array_16, s.nptr);
                        m_s16.scan(m_edges);
                        break;
                    }
                    case MART_NODE_32: {
                        m_s32.set(&m_array_32, s.nptr);
                        m_s32.scan(m_edges);
                        break;
                    }
                    case MART_NODE_64: {
                        m_s64.set(&m_array_64, s.nptr);
                        for (int k = 0; k <= r; ++k) {
                            for (int i = bpt[k]; i < bpt[k + 1]; ++i) {
                                DEBUG_ABORT_IF_NE(k, hdt[lut[i]]);
                                const mart_pointer cptr = m_s64.find(lut[i]);
                                if (!check_mart_nullptr(cptr)) {
                                    m_edges.push_back(mart_edge{lut[i], cptr});
                                }
                            }
                        }
                        break;
                    }
                    case MART_NODE_128: {
                        m_s128.set(&m_array_128, s.nptr);
                        for (int k = 0; k <= r; ++k) {
                            for (int i = bpt[k]; i < bpt[k + 1]; ++i) {
                                DEBUG_ABORT_IF_NE(k, hdt[lut[i]]);
                                const mart_pointer cptr = m_s128.find(lut[i]);
                                if (!check_mart_nullptr(cptr)) {
                                    m_edges.push_back(mart_edge{lut[i], cptr});
                                }
                            }
                        }
                        break;
                    }
                    case MART_NODE_256: {
                        m_s256.set(&m_array_256, s.nptr);
                        for (int k = 0; k <= r; ++k) {
                            for (int i = bpt[k]; i < bpt[k + 1]; ++i) {
                                DEBUG_ABORT_IF_NE(k, hdt[lut[i]]);
                                const mart_pointer cptr = m_s256.find(lut[i]);
                                if (!check_mart_nullptr(cptr)) {
                                    m_edges.push_back(mart_edge{lut[i], cptr});
                                }
                            }
                        }
                        break;
                    }
                    default: {
                        ABORT_IF(true);  // should not come
                        break;
                    }
                }
#else
                if (s.nptr.ntype == MART_NODE_4) {
                    m_s4.set(&m_array_4, s.nptr);
                    m_s4.scan(m_edges);
                } else if (s.nptr.ntype == MART_NODE_32) {
                    m_s32.set(&m_array_32, s.nptr);
                    m_s32.scan(m_edges);
                } else if (s.nptr.ntype == MART_NODE_256) {
                    m_s256.set(&m_array_256, s.nptr);
                    for (int k = 0; k <= r; ++k) {
                        for (int i = bpt[k]; i < bpt[k + 1]; ++i) {
                            DEBUG_ABORT_IF_NE(k, hdt[lut[i]]);
                            const mart_pointer cptr = m_s256.find(lut[i]);
                            if (!check_mart_nullptr(cptr)) {
                                m_edges.push_back(mart_edge{lut[i], cptr});
                            }
                        }
                    }
                } else {
                    ABORT_IF(true);  // should not come
                }
#endif
                for (const mart_edge& e : m_edges) {
                    if (hdt[e.c] <= r) {
                        if (e.nptr.ntype == MART_LEAF) {  // leaf
                            const uint32_t lpos = e.nptr.nid;
                            auto [bptr, eptr] = m_idseqs.access(lpos);
                            for (auto it = bptr; it != eptr; ++it) {
                                fn(*it);
                            }
                        } else {  // internal
                            m_states.push_back(state_type{e.nptr, s.bpos + 1, s.dist + hdt[e.c]});
                        }
                    }
                }
            } else {  // exact match
                while (true) {
                    const uint8_t c = m_chunks[s.bpos];

#ifdef MART_SPACE_EFFICIENT
                    switch (s.nptr.ntype) {
                        case MART_NODE_2: {
                            s.nptr = m_array_2.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_4: {
                            s.nptr = m_array_4.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_8: {
                            s.nptr = m_array_8.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_16: {
                            s.nptr = m_array_16.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_32: {
                            s.nptr = m_array_32.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_64: {
                            s.nptr = m_array_64.find_child(s.nptr, c);
                            break;
                        }
                        case MART_NODE_128: {
                            s.nptr = m_array_128.find_child(s.nptr, c);
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
#else
                    if (s.nptr.ntype == MART_NODE_4) {
                        s.nptr = m_array_4.find_child(s.nptr, c);
                    } else if (s.nptr.ntype == MART_NODE_32) {
                        s.nptr = m_array_32.find_child(s.nptr, c);
                    } else if (s.nptr.ntype == MART_NODE_256) {
                        s.nptr = m_array_256.find_child(s.nptr, c);
                    } else {
                        ABORT_IF(true);  // should not come
                    }
#endif

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

                    s.bpos += 1;
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
        return m_infacts_ptr and (get_ls_cost() < get_trie_cost());
    }

    size_t get_split_count() override {
        return m_split_count;
    }

    void innode_stats() override {
        tfm::printfln("## innode_stats ##");
#ifdef MART_SPACE_EFFICIENT
        tfm::printfln("- num node 2:   %d (%d)", m_array_2.num_nodes(), m_array_2.num_emp_nodes());
        tfm::printfln("- num node 4:   %d (%d)", m_array_4.num_nodes(), m_array_4.num_emp_nodes());
        tfm::printfln("- num node 8:   %d (%d)", m_array_8.num_nodes(), m_array_8.num_emp_nodes());
        tfm::printfln("- num node 16:  %d (%d)", m_array_16.num_nodes(), m_array_16.num_emp_nodes());
        tfm::printfln("- num node 32:  %d (%d)", m_array_32.num_nodes(), m_array_32.num_emp_nodes());
        tfm::printfln("- num node 64:  %d (%d)", m_array_64.num_nodes(), m_array_64.num_emp_nodes());
        tfm::printfln("- num node 128: %d (%d)", m_array_128.num_nodes(), m_array_128.num_emp_nodes());
        tfm::printfln("- num node 256: %d (%d)", m_array_256.num_nodes(), m_array_256.num_emp_nodes());
#else
        tfm::printfln("- num node 4:   %d (%d)", m_array_4.num_nodes(), m_array_4.num_emp_nodes());
        tfm::printfln("- num node 32:  %d (%d)", m_array_32.num_nodes(), m_array_32.num_emp_nodes());
        tfm::printfln("- num node 256: %d (%d)", m_array_256.num_nodes(), m_array_256.num_emp_nodes());
#endif
    }

    void population_stats() override {
#ifdef MART_SPACE_EFFICIENT
        m_array_2.population_stats();
        m_array_4.population_stats();
        m_array_8.population_stats();
        m_array_16.population_stats();
        m_array_32.population_stats();
        m_array_64.population_stats();
        m_array_128.population_stats();
        // m_array_256.population_stats();
#else
        m_array_4.population_stats();
        m_array_32.population_stats();
#endif
    }

    void debug_arrays() {
#ifdef MART_SPACE_EFFICIENT
#else
        m_array_4.debug_emplist();
        m_array_32.debug_emplist();
#endif
    }

  private:
    mart_insert_flags expand(mart_cursor& mc, uint8_t c, mart_pointer new_ptr) {
        m_edges.clear();

#ifdef MART_SPACE_EFFICIENT
        switch (mc.nptr.ntype) {
            case MART_NODE_2: {
                m_array_2.extract(mc, m_edges);
                mc.nptr = m_array_4.make_node(m_edges);
                update_srcptr(mc);
                return m_array_4.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_4: {
                m_array_4.extract(mc, m_edges);
                mc.nptr = m_array_8.make_node(m_edges);
                update_srcptr(mc);
                return m_array_8.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_8: {
                m_array_8.extract(mc, m_edges);
                mc.nptr = m_array_16.make_node(m_edges);
                update_srcptr(mc);
                return m_array_16.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_16: {
                m_array_16.extract(mc, m_edges);
                mc.nptr = m_array_32.make_node(m_edges);
                update_srcptr(mc);
                return m_array_32.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_32: {
                m_array_32.extract(mc, m_edges);
                mc.nptr = m_array_64.make_node(m_edges);
                update_srcptr(mc);
                return m_array_64.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_64: {
                m_array_64.extract(mc, m_edges);
                mc.nptr = m_array_128.make_node(m_edges);
                update_srcptr(mc);
                return m_array_128.append_ptr(mc, c, new_ptr);
            }
            case MART_NODE_128: {
                m_array_128.extract(mc, m_edges);
                mc.nptr = m_array_256.make_node(m_edges);
                update_srcptr(mc);
                return m_array_256.insert_ptr(mc, c, new_ptr);
            }
            default: {
                ABORT_IF(true);  // should not come
                break;
            }
        }
#else
        if (mc.nptr.ntype == MART_NODE_4) {
            m_array_4.extract(mc, m_edges);
            mc.nptr = m_array_32.make_node(m_edges);
            update_srcptr(mc);
            return m_array_32.append_ptr(mc, c, new_ptr);
        } else if (mc.nptr.ntype == MART_NODE_32) {
            m_array_32.extract(mc, m_edges);
            mc.nptr = m_array_256.make_node(m_edges);
            update_srcptr(mc);
            return m_array_256.insert_ptr(mc, c, new_ptr);
        } else {
            ABORT_IF(true);  // should not come
        }
#endif
    }

    void update_srcptr(const mart_cursor& mc) {
#ifdef MART_SPACE_EFFICIENT
        switch (mc.pptr.ntype) {
            case MART_NODE_2: {
                m_array_2.update_srcptr(mc);
                break;
            }
            case MART_NODE_4: {
                m_array_4.update_srcptr(mc);
                break;
            }
            case MART_NODE_8: {
                m_array_8.update_srcptr(mc);
                break;
            }
            case MART_NODE_16: {
                m_array_16.update_srcptr(mc);
                break;
            }
            case MART_NODE_32: {
                m_array_32.update_srcptr(mc);
                break;
            }
            case MART_NODE_64: {
                m_array_64.update_srcptr(mc);
                break;
            }
            case MART_NODE_128: {
                m_array_128.update_srcptr(mc);
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
#else
        if (mc.pptr.ntype == MART_NODE_4) {
            m_array_4.update_srcptr(mc);
        } else if (mc.pptr.ntype == MART_NODE_32) {
            m_array_32.update_srcptr(mc);
        } else if (mc.pptr.ntype == MART_NODE_256) {
            m_array_256.update_srcptr(mc);
        } else {
            ABORT_IF(true);  // should not come
        }
#endif
    }

    void split_node(mart_cursor& mc, const int bpos) {
        ABORT_IF_LE(bpos, m_bpos_beg);
        ABORT_IF_NE(mc.nptr.ntype, MART_LEAF);

        m_split_count += 1;

        const uint32_t lpos = mc.nptr.nid;
        ABORT_IF_LE(m_idseqs.size(), lpos);

        // update cost
        if (m_infacts_ptr) {
            const int _depth = (bpos - m_bpos_beg - 1) * m_ints_per_chunk;
            const int g_size = m_idseqs.group_size(lpos);
            m_trie_cost -= g_size * m_outfacts_ptr[_depth];
            m_trie_cost += (m_infacts_ptr[_depth] * m_in_weight) + (g_size * m_outfacts_ptr[_depth + m_ints_per_chunk]);
        }

        std::vector<uint32_t> idbufs[256];
        {
            const auto ids = m_idseqs.extract(lpos);
            for (uint32_t id : ids) {
                const vint_type* vcode = m_database->access(id);
                idbufs[to_chunk(vcode, bpos)].push_back(id);
            }
        }

        m_edges.clear();

        for (int c = 0; c < 256; c++) {
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

#ifdef MART_SPACE_EFFICIENT
        if (m_edges.size() <= 2) {
            mc.nptr = m_array_2.make_node(m_edges);
        } else if (m_edges.size() <= 4) {
            mc.nptr = m_array_4.make_node(m_edges);
        } else if (m_edges.size() <= 8) {
            mc.nptr = m_array_8.make_node(m_edges);
        } else if (m_edges.size() <= 16) {
            mc.nptr = m_array_16.make_node(m_edges);
        } else if (m_edges.size() <= 32) {
            mc.nptr = m_array_32.make_node(m_edges);
        } else if (m_edges.size() <= 64) {
            mc.nptr = m_array_64.make_node(m_edges);
        } else if (m_edges.size() <= 128) {
            mc.nptr = m_array_128.make_node(m_edges);
        } else {
            mc.nptr = m_array_256.make_node(m_edges);
        }
#else
        if (m_edges.size() <= 4) {
            mc.nptr = m_array_4.make_node(m_edges);
        } else if (m_edges.size() <= 32) {
            mc.nptr = m_array_32.make_node(m_edges);
        } else {
            mc.nptr = m_array_256.make_node(m_edges);
        }
#endif

        update_srcptr(mc);
    }

    uint8_t to_chunk(const vint_type* vcode, int bpos) {
        if (get_bits() == 1) {
            return *vcode >> (bpos * m_ints_per_chunk);
        }

        const int s = bpos * m_ints_per_chunk;
        const int e = std::min((bpos + 1) * m_ints_per_chunk, N);

        uint8_t chunk = 0;
        for (int i = s; i < e; i++) {
            const int shift = (i - s) * get_bits();
            chunk |= vcode_tools<N>::get_int(vcode, i, get_bits()) << shift;
        }
        return chunk;
    }

    void to_chunks(const vint_type* vcode) {
        for (int bpos = m_bpos_beg; bpos < m_bpos_end; bpos++) {
            m_chunks[bpos] = to_chunk(vcode, bpos);
        }
    }
};

}  // namespace dyft
