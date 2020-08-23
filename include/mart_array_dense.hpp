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
class mart_array_dense {
    static_assert(2 <= K and K < 255);

  public:
    static constexpr uint64_t PTR_SIZE = 5;
    static constexpr uint64_t IDXS_OFFSET = 1;
    static constexpr uint64_t PTRS_OFFSET = 1 + 256;

    static constexpr uint64_t BYTES = 1 + 256 + (K * PTR_SIZE);  // header + 256 indexes + K ptrs
    static constexpr mart_node_types NTYPE = NodeType;

    static constexpr uint8_t NIL_IDX = UINT8_MAX;

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
    mart_array_dense() = default;

    mart_insert_flags insert_ptr(mart_cursor& mc, uint8_t c, mart_pointer new_ptr) {
        DEBUG_ABORT_IF_NE(mc.nptr.ntype, NTYPE);

        const uint64_t pos = mc.nptr.nid * BYTES;

        const uint32_t num = m_nodes[pos];
        uint8_t* idxs = &m_nodes[pos + IDXS_OFFSET];
        uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        DEBUG_ABORT_IF_LT(K, num);

        if (idxs[c] != NIL_IDX) {  // found?
            const uint32_t i = idxs[c];
            mc = mart_cursor{i, mc.nptr, get_martptr(ptrs + (i * PTR_SIZE))};
            return MART_FOUND;
        }

        // can insert
        if (num < K) {
            m_nodes[pos] += 1;
            idxs[c] = static_cast<uint8_t>(num);
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
        uint8_t* idxs = &m_nodes[pos + IDXS_OFFSET];
        uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        DEBUG_ABORT_IF_LE(K, num);
        DEBUG_ABORT_IF_NE(idxs[c], NIL_IDX);

        // not searched
        m_nodes[pos] += 1;
        idxs[c] = static_cast<uint8_t>(num);
        set_martptr(ptrs + (num * PTR_SIZE), new_ptr);

        // set src
        mc = mart_cursor{num, mc.nptr, new_ptr};
        return MART_INSERTED;
    }

    void extract(const mart_cursor& mc, std::vector<mart_edge>& edges) {
        DEBUG_ABORT_IF_NE(mc.nptr.ntype, NTYPE);

        const uint32_t nid = mc.nptr.nid;
        const uint64_t pos = nid * BYTES;

        const uint8_t* idxs = &m_nodes[pos + IDXS_OFFSET];
        const uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        for (uint32_t i = 0; i < 256; i++) {
            if (idxs[i] != NIL_IDX) {
                edges.push_back(mart_edge{uint8_t(i), get_martptr(ptrs + (idxs[i] * PTR_SIZE))});
            }
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

            const uint64_t new_pos = new_nid * BYTES;

            m_nodes[new_pos] = 0;  // fanout
            uint8_t* idxs = &m_nodes[new_pos + IDXS_OFFSET];
            std::fill(idxs, idxs + 256, NIL_IDX);

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

        uint8_t* idxs = &m_nodes[new_pos + IDXS_OFFSET];
        std::fill(idxs, idxs + 256, NIL_IDX);

        return mart_pointer{new_nid, NTYPE};
    }

    // returns new nid
    mart_pointer make_node(const std::vector<mart_edge>& edges) {
        ABORT_IF_OUT(edges.size(), 0, K);

        const mart_pointer new_nptr = make_node();

        const uint32_t new_nid = new_nptr.nid;
        const uint64_t new_pos = new_nid * BYTES;

        uint8_t* idxs = &m_nodes[new_pos + IDXS_OFFSET];
        uint8_t* ptrs = &m_nodes[new_pos + PTRS_OFFSET];

        m_nodes[new_pos] = static_cast<uint8_t>(edges.size());
        for (uint32_t i = 0; i < edges.size(); i++) {
            idxs[edges[i].c] = static_cast<uint8_t>(i);
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

        const uint8_t* idxs = &m_nodes[pos + IDXS_OFFSET];
        const uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        if (idxs[c] != NIL_IDX) {  // found?
            return get_martptr(ptrs + (idxs[c] * PTR_SIZE));
        }
        return make_mart_nullptr();
    }

    struct searcher {
        searcher() = default;

        searcher(const mart_array_dense* obj, mart_pointer nptr) {
            set(obj, nptr);
        }

        void set(const mart_array_dense* obj, mart_pointer nptr) {
            DEBUG_ABORT_IF_NE(nptr.ntype, NTYPE);

            const uint64_t pos = nptr.nid * BYTES;

            m_num = obj->m_nodes[pos];
            m_idxs = &obj->m_nodes[pos + IDXS_OFFSET];
            m_ptrs = &obj->m_nodes[pos + PTRS_OFFSET];
        }

        mart_pointer find(uint8_t c) {
            if (m_idxs[c] != NIL_IDX) {  // found?
                return get_martptr(m_ptrs + (m_idxs[c] * PTR_SIZE));
            }
            return make_mart_nullptr();
        }
        void scan(std::vector<mart_edge>& edges) {
            for (uint32_t i = 0; i < 256; i++) {
                if (m_idxs[i] != NIL_IDX) {
                    edges.push_back(mart_edge{uint8_t(i), get_martptr(m_ptrs + (m_idxs[i] * PTR_SIZE))});
                }
            }
        }

      private:
        uint32_t m_num = 0;
        const uint8_t* m_idxs = nullptr;
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

        tfm::printfln("## mart_array_dense (K=%d) ##", K);
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
};

}  // namespace dyft
