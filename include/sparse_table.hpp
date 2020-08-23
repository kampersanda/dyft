#pragma once

#include "sparse_group.hpp"

namespace dyft {

class sparse_table {
  private:
    std::vector<sparse_group> m_groups;
    uint32_t m_size = 0;

  public:
    sparse_table() = default;

    std::pair<uint32_t*, uint32_t*> access(uint32_t idx) {
        ABORT_IF_LE(m_size, idx);
        const uint32_t gpos = idx / sparse_group::SIZE;
        const uint32_t gmod = idx % sparse_group::SIZE;
        return m_groups[gpos].access(gmod);
    }

    void push_back() {
        if (m_size / sparse_group::SIZE == m_groups.size()) {
            m_groups.push_back(sparse_group());
        }
        m_size += 1;
    }

    void push_back(const std::vector<uint32_t>& datvec) {
        if (m_size / sparse_group::SIZE == m_groups.size()) {
            m_groups.push_back(sparse_group());
        }
        const uint32_t gpos = m_size / sparse_group::SIZE;
        const uint32_t gmod = m_size % sparse_group::SIZE;
        m_groups[gpos].insert(gmod, datvec);
        m_size += 1;
    }

    void insert(uint32_t idx, uint32_t data) {
        ABORT_IF_LE(m_size, idx);
        const uint32_t gpos = idx / sparse_group::SIZE;
        const uint32_t gmod = idx % sparse_group::SIZE;
        m_groups[gpos].insert(gmod, data);
    }

    void insert(uint32_t idx, const std::vector<uint32_t>& datvec) {
        ABORT_IF_LE(m_size, idx);
        const uint32_t gpos = idx / sparse_group::SIZE;
        const uint32_t gmod = idx % sparse_group::SIZE;
        m_groups[gpos].insert(gmod, datvec);
    }

    std::vector<uint32_t> extract(uint32_t idx) {
        ABORT_IF_LE(m_size, idx);
        const uint32_t gpos = idx / sparse_group::SIZE;
        const uint32_t gmod = idx % sparse_group::SIZE;
        return m_groups[gpos].extract(gmod);
    }

    uint32_t size() const {
        return m_size;
    }

    uint32_t group_size(uint32_t idx) const {
        ABORT_IF_LE(m_size, idx);
        const uint32_t gpos = idx / sparse_group::SIZE;
        const uint32_t gmod = idx % sparse_group::SIZE;
        return m_groups[gpos].size(gmod);
    }
};

}  // namespace dyft
