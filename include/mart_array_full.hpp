#pragma once

#include <vector>

#include "abort_if.hpp"
#include "mart_common.hpp"

namespace dyft {

template <mart_node_types NodeType>
class mart_array_full {
  public:
    static constexpr uint64_t PTR_SIZE = 5;
    static constexpr uint64_t PTRS_OFFSET = 1;

    static constexpr uint64_t BYTES = 1 + (256 * PTR_SIZE);  // header + 256 ptrs
    static constexpr mart_node_types NTYPE = NodeType;

  private:
    std::vector<uint8_t> m_nodes;

    static mart_pointer get_martptr(const uint8_t* ptrs) {
        return {*reinterpret_cast<const uint32_t*>(ptrs), mart_node_types(ptrs[4])};
    }
    static void set_martptr(uint8_t* ptrs, mart_pointer new_ptr) {
        *reinterpret_cast<uint32_t*>(ptrs) = new_ptr.nid;
        ptrs[4] = uint8_t(new_ptr.ntype);
    }

  public:
    mart_array_full() = default;

    mart_insert_flags insert_ptr(mart_cursor& mc, uint8_t c, mart_pointer new_ptr) {
        DEBUG_ABORT_IF_NE(mc.nptr.ntype, NTYPE);

        const uint64_t pos = mc.nptr.nid * BYTES;

        uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];
        mart_pointer mptr = get_martptr(ptrs + (c * PTR_SIZE));

        if (mptr.nid == MART_NILID) {
            // not found
            m_nodes[pos] += 1;
            set_martptr(ptrs + (c * PTR_SIZE), new_ptr);
            mc = mart_cursor{uint32_t(c), mc.nptr, new_ptr};
            return MART_INSERTED;
        }

        // found
        mc = mart_cursor{uint32_t(c), mc.nptr, mptr};
        return MART_FOUND;
    }

    // returns new node ptr
    mart_pointer make_node() {
        if (m_nodes.size() == 0) {
            m_nodes.reserve(BYTES);
        } else if (m_nodes.size() + BYTES > m_nodes.capacity()) {
            m_nodes.reserve(m_nodes.capacity() * 2);
        }

        const uint64_t new_pos = m_nodes.size();
        const uint32_t new_nid = static_cast<uint32_t>(new_pos / BYTES);
        m_nodes.resize(new_pos + BYTES);
        m_nodes[new_pos] = 0;  // fanout

        uint8_t* ptrs = &m_nodes[new_pos + PTRS_OFFSET];
        for (uint32_t i = 0; i < 256; i++) {
            set_martptr(ptrs + (i * PTR_SIZE), make_mart_nullptr());
        }

        return mart_pointer{new_nid, NTYPE};
    }

    // returns new nid
    mart_pointer make_node(const std::vector<mart_edge>& edges) {
        const mart_pointer new_nptr = make_node();

        const uint32_t new_nid = new_nptr.nid;
        const uint64_t new_pos = new_nid * BYTES;

        m_nodes[new_pos] = static_cast<uint8_t>(edges.size());
        uint8_t* ptrs = &m_nodes[new_pos + PTRS_OFFSET];

        for (const mart_edge& e : edges) {
            set_martptr(ptrs + (e.c * PTR_SIZE), e.nptr);
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

        const uint64_t pos = nptr.nid * BYTES;
        const uint8_t* ptrs = &m_nodes[pos + PTRS_OFFSET];

        return get_martptr(ptrs + (c * PTR_SIZE));
    }

    struct searcher {
        searcher() = default;

        searcher(const mart_array_full* obj, mart_pointer nptr) {
            set(obj, nptr);
        }

        void set(const mart_array_full* obj, mart_pointer nptr) {
            DEBUG_ABORT_IF_NE(nptr.ntype, NTYPE);

            const uint64_t pos = nptr.nid * BYTES;
            m_ptrs = &obj->m_nodes[pos + PTRS_OFFSET];
        }

        mart_pointer find(uint8_t c) {
            return get_martptr(m_ptrs + (c * PTR_SIZE));
        }
        void scan(std::vector<mart_edge>& edges) {
            for (uint32_t i = 0; i < 256; i++) {
                mart_pointer ptr = get_martptr(m_ptrs + (i * PTR_SIZE));
                if (ptr.nid != MART_NILID) {
                    edges.push_back(mart_edge{uint8_t(i), ptr});
                }
            }
        }

      private:
        const uint8_t* m_ptrs = nullptr;
    };

    uint32_t num_nodes() {
        return m_nodes.size() / BYTES;
    }
    uint32_t num_emp_nodes() {
        return 0;
    }
};

}  // namespace dyft
