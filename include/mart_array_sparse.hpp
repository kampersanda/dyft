#pragma once

#if defined(__AVX__) || defined(__AVX2__)
#include <immintrin.h>
#endif

#include <set>
#include <vector>

#include "abort_if.hpp"
#include "mart_common.hpp"

namespace dyft {

template <uint32_t K, mart_node_types NodeType>
class mart_array_sparse {
    static_assert(2 <= K and K <= 255);

  public:
    static constexpr uint64_t PTR_SIZE = 5;
    static constexpr uint64_t KEYS_OFFSET = 1;
    static constexpr uint64_t PTRS_OFFSET = 1 + K;

    static constexpr uint64_t BYTES = 1 + K + (K * PTR_SIZE);  // header + K keys + K ptrs
    static constexpr mart_node_types NTYPE = NodeType;

  private:
    std::vector<uint8_t> m_nodes;
    uint32_t m_head_nid = MART_NILID;
    uint32_t m_num_emps = 0;

    static mart_pointer get_martptr(const uint8_t* ptrs) {
        return {*reinterpret_cast<const uint32_t*>(ptrs), mart_node_types(ptrs[4])};
    }
    static void set_martptr(uint8_t* ptrs, mart_pointer new_ptr) {
        *reinterpret_cast<uint32_t*>(ptrs) = new_ptr.nid;
        ptrs[4] = uint8_t(new_ptr.ntype);
    }

  public:
    mart_array_sparse() = default;

    mart_insert_flags insert_ptr(mart_cursor& mc, uint8_t c, mart_pointer new_ptr) {
        DEBUG_ABORT_IF_NE(mc.nptr.ntype, NTYPE);

        const uint64_t pos = mc.nptr.nid * BYTES;

        const uint32_t num = m_nodes[pos];
        uint8_t* keys = &m_nodes[pos + KEYS_OFFSET];
        uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        DEBUG_ABORT_IF_LT(K, num);

        for (uint32_t i = 0; i < num; i++) {
            if (keys[i] == c) {  // found?
                mc = mart_cursor{i, mc.nptr, get_martptr(ptrs + (i * PTR_SIZE))};
                return MART_FOUND;
            }
        }

        // can insert
        if (num < K) {
            m_nodes[pos] += 1;
            keys[num] = c;
            set_martptr(ptrs + (num * PTR_SIZE), new_ptr);

            mc = mart_cursor{num, mc.nptr, new_ptr};
            return MART_INSERTED;  // inserted
        }

        return MART_NEEDED_TO_EXPAND;
    }

    mart_insert_flags append_ptr(mart_cursor& mc, uint8_t c, mart_pointer new_ptr) {
        DEBUG_ABORT_IF_NE(mc.nptr.ntype, NTYPE);

        const uint64_t pos = mc.nptr.nid * BYTES;

        const uint32_t num = m_nodes[pos];
        uint8_t* keys = &m_nodes[pos + KEYS_OFFSET];
        uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        DEBUG_ABORT_IF_LE(K, num);

        // not searched
        m_nodes[pos] += 1;
        keys[num] = c;
        set_martptr(ptrs + (num * PTR_SIZE), new_ptr);

        // set src
        mc = mart_cursor{num, mc.nptr, new_ptr};
        return MART_INSERTED;
    }

    void extract(const mart_cursor& mc, std::vector<mart_edge>& edges) {
        DEBUG_ABORT_IF_NE(mc.nptr.ntype, NTYPE);

        const uint32_t nid = mc.nptr.nid;
        const uint64_t pos = nid * BYTES;

        const uint32_t num = m_nodes[pos];
        const uint8_t* keys = &m_nodes[pos + KEYS_OFFSET];
        const uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        for (uint32_t i = 0; i < num; i++) {
            edges.push_back(mart_edge{keys[i], get_martptr(ptrs + (i * PTR_SIZE))});
        }

        // To make empty
        if (m_head_nid == MART_NILID) {
            prev_ref(nid) = nid;
            next_ref(nid) = nid;
            m_head_nid = nid;
        } else {
            const uint32_t prev_nid = prev_ref(m_head_nid);
            prev_ref(nid) = prev_nid;
            next_ref(nid) = m_head_nid;
            next_ref(prev_nid) = nid;
            prev_ref(m_head_nid) = nid;
        }

        m_num_emps += 1;
    }

    mart_pointer make_node() {
        // Reuse empty element
        if (m_head_nid != MART_NILID) {
            m_num_emps -= 1;

            const uint32_t new_nid = m_head_nid;
            const uint32_t prev_nid = prev_ref(m_head_nid);
            const uint32_t next_nid = next_ref(m_head_nid);

            if (next_nid == m_head_nid) {
                m_head_nid = MART_NILID;
            } else {
                m_head_nid = next_nid;
                prev_ref(next_nid) = prev_nid;
                next_ref(prev_nid) = next_nid;
            }
            m_nodes[new_nid * BYTES] = 0;  // fanout
            return mart_pointer{new_nid, NTYPE};
        }

        if (m_nodes.size() == 0) {
            m_nodes.reserve(BYTES);
        } else if (m_nodes.size() + BYTES > m_nodes.capacity()) {
            m_nodes.reserve(m_nodes.capacity() * 2);
        }

        const uint64_t new_pos = m_nodes.size();
        const uint32_t new_nid = static_cast<uint32_t>(new_pos / BYTES);
        m_nodes.resize(new_pos + BYTES);
        m_nodes[new_pos] = 0;  // fanout

        return mart_pointer{new_nid, NTYPE};
    }

    // returns new nid
    mart_pointer make_node(const std::vector<mart_edge>& edges) {
        ABORT_IF_OUT(edges.size(), 0, K);

        const mart_pointer new_nptr = make_node();

        const uint32_t new_nid = new_nptr.nid;
        const uint64_t new_pos = new_nid * BYTES;

        uint8_t* keys = &m_nodes[new_pos + KEYS_OFFSET];
        uint8_t* ptrs = &m_nodes[new_pos + PTRS_OFFSET];

        m_nodes[new_pos] = static_cast<uint8_t>(edges.size());
        for (uint32_t i = 0; i < edges.size(); i++) {
            keys[i] = edges[i].c;
            set_martptr(ptrs + (i * PTR_SIZE), edges[i].nptr);
        }

        return new_nptr;
    }

    void update_srcptr(const mart_cursor& mc) {
        DEBUG_ABORT_IF_NE(mc.pptr.ntype, NTYPE);

        const uint64_t pos = mc.pptr.nid * BYTES;
        uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        set_martptr(ptrs + (mc.offset * PTR_SIZE), mc.nptr);
    }

    mart_pointer find_child(mart_pointer nptr, uint8_t c) {
        DEBUG_ABORT_IF_NE(nptr.ntype, NTYPE);

        const uint32_t nid = nptr.nid;
        const uint64_t pos = nid * BYTES;

        const uint32_t num = m_nodes[pos];
        const uint8_t* keys = &m_nodes[pos + KEYS_OFFSET];
        const uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        return adaptive_find(num, keys, ptrs, c);
    }

    struct searcher {
        searcher() = default;

        searcher(const mart_array_sparse* obj, mart_pointer nptr) {
            set(obj, nptr);
        }

        void set(const mart_array_sparse* obj, mart_pointer nptr) {
            DEBUG_ABORT_IF_NE(nptr.ntype, NTYPE);

            const uint64_t pos = nptr.nid * BYTES;

            m_num = obj->m_nodes[pos];
            m_keys = &obj->m_nodes[pos + KEYS_OFFSET];
            m_ptrs = &obj->m_nodes[pos + PTRS_OFFSET];
        }

        mart_pointer find(uint8_t c) {
            return adaptive_find(m_num, m_keys, m_ptrs, c);
        }
        void scan(std::vector<mart_edge>& edges) {
            for (uint32_t i = 0; i < m_num; i++) {
                edges.push_back(mart_edge{m_keys[i], get_martptr(m_ptrs + (i * PTR_SIZE))});
            }
        }

      private:
        uint32_t m_num = 0;
        const uint8_t* m_keys = nullptr;
        const uint8_t* m_ptrs = nullptr;
    };

    void population_stats() {
        std::set<uint64_t> emps = get_emps();

        uint64_t cnts[K + 1];
        std::fill(cnts, cnts + K + 1, 0);

        for (uint64_t i = 0; i < m_nodes.size(); i += BYTES) {
            if (emps.find(i) != emps.end()) {
                continue;
            }
            DEBUG_ABORT_IF_LT(K, m_nodes[i]);
            cnts[m_nodes[i]] += 1;
        }

        const double sum = std::accumulate(cnts, cnts + K + 1, 0.0);

        tfm::printfln("## mart_array_sparse (K=%d) ##", K);
        for (uint32_t k = 0; k <= K; k++) {
            if (cnts[k] == 0) continue;
            tfm::printfln("- k=%d: %d (%g)", k, cnts[k], cnts[k] / sum);
        }
    }

    uint32_t num_nodes() {
        return (m_nodes.size() / BYTES) - m_num_emps;
    }
    uint32_t num_emp_nodes() {
        return m_num_emps;
    }

    void debug_emplist() {
        uint32_t emps = 0;
        if (m_head_nid != MART_NILID) {
            uint32_t emp_nid = m_head_nid;
            while (true) {
                emps += 1;
                uint32_t next_nid = next_ref(emp_nid);
                if (next_nid == m_head_nid) break;
                emp_nid = next_nid;
            }
        }
        ABORT_IF_NE(m_num_emps, emps);
    }

    std::set<uint64_t> get_emps() {
        std::set<uint64_t> emps;
        if (m_head_nid != MART_NILID) {
            uint32_t emp_nid = m_head_nid;
            while (true) {
                emps.insert(emp_nid * BYTES);
                uint32_t next_nid = next_ref(emp_nid);
                if (next_nid == m_head_nid) break;
                emp_nid = next_nid;
            }
        }
        return emps;
    }

  private:
    uint32_t& prev_ref(uint32_t nid) {
        return reinterpret_cast<uint32_t*>(&m_nodes[nid * BYTES])[0];
    }
    uint32_t& next_ref(uint32_t nid) {
        return reinterpret_cast<uint32_t*>(&m_nodes[nid * BYTES])[1];
    }

    static mart_pointer adaptive_find(uint32_t num, const uint8_t* keys, const uint8_t* ptrs, uint8_t c) {
        if constexpr (K <= 8) {
            // Linear Scan
            for (uint32_t i = 0; i < num; i++) {
                if (keys[i] == c) {  // found?
                    return get_martptr(ptrs + (i * PTR_SIZE));
                }
            }
        } else {
#ifdef __AVX2__
            __m256i c_256i = _mm256_set1_epi8(c);
            for (uint32_t i = 0; i < num; i += 32) {
                __m256i cmp = _mm256_cmpeq_epi8(c_256i, _mm256_loadu_si256((__m256i*)(keys + i)));
                const unsigned check_bits = _mm256_movemask_epi8(cmp) & ((1ULL << std::min(32U, num - i)) - 1ULL);
                if (check_bits) {
                    const unsigned j = __builtin_ctz(check_bits) + i;
                    return get_martptr(ptrs + (j * PTR_SIZE));
                }
            }
#elif __AVX__
            __m128i c_128i = _mm_set1_epi8(c);
            for (uint32_t i = 0; i < num; i += 16) {
                __m128i cmp = _mm_cmpeq_epi8(c_128i, _mm_loadu_si128((__m128i*)(keys + i)));
                const unsigned check_bits = _mm_movemask_epi8(cmp) & ((1U << std::min(16U, num - i)) - 1U);
                if (check_bits) {
                    const unsigned j = __builtin_ctz(check_bits) + i;
                    return get_martptr(ptrs + (j * PTR_SIZE));
                }
            }
#else
            // Linear Scan
            for (uint32_t i = 0; i < num; i++) {
                if (keys[i] == c) {  // found?
                    return get_martptr(ptrs + (i * PTR_SIZE));
                }
            }
#endif
        }
        return make_mart_nullptr();
    }
};

}  // namespace dyft
