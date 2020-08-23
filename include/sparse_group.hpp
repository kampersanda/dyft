#pragma once

#include "abort_if.hpp"
#include "bit_tools.hpp"

namespace dyft {

class sparse_group {
  public:
    static constexpr uint32_t SIZE = 64;

  private:
    uint64_t m_bitmap = 0;
    std::vector<uint32_t> m_group;

  public:
    sparse_group() = default;

    std::pair<uint32_t*, uint32_t*> access(uint32_t idx) {
        ABORT_IF_LE(SIZE, idx);

        if ((m_bitmap & (1ULL << idx)) == 0ULL) {
            return {nullptr, nullptr};
        }

        const uint64_t bitmask = (1ULL << idx) - 1ULL;
        const int howmany = bit_tools::popcnt(m_bitmap & bitmask);
        const int totones = bit_tools::popcnt(m_bitmap);

        const uint32_t size = m_group[howmany + 1] - m_group[howmany];
        uint32_t* ptr = m_group.data() + (totones + 1 + m_group[howmany]);

        return {ptr, ptr + size};
    }

    void insert(uint32_t idx, uint32_t data) {
        ABORT_IF_LE(SIZE, idx);

        if (m_bitmap == 0) {
            m_bitmap = 1ULL << idx;
            m_group = std::vector<uint32_t>{0, 1, data};
            return;
        }

        const uint64_t bitmask = (1ULL << idx) - 1ULL;
        const int howmany = bit_tools::popcnt(m_bitmap & bitmask);

        if ((m_bitmap & (1ULL << idx)) == 0ULL) {
            m_group.insert(m_group.begin() + howmany, m_group[howmany]);
            m_bitmap |= (1ULL << idx);
        }

        const int totones = bit_tools::popcnt(m_bitmap);

        // totones + 1 = dataの始点
        // m_group[howmany + 1] = idxの末尾
        m_group.insert(m_group.begin() + (totones + 1 + m_group[howmany + 1]), data);
        for (int i = howmany + 1; i <= totones; i++) {
            m_group[i] += 1;
        }
    }

    void insert(uint32_t idx, const std::vector<uint32_t>& datvec) {
        ABORT_IF_LE(SIZE, idx);

        if (m_bitmap == 0) {
            m_bitmap = 1ULL << idx;
            m_group = std::vector<uint32_t>{0, static_cast<uint32_t>(datvec.size())};
            std::copy(datvec.begin(), datvec.end(), std::back_inserter(m_group));
            return;
        }

        const uint64_t bitmask = (1ULL << idx) - 1ULL;
        const int howmany = bit_tools::popcnt(m_bitmap & bitmask);

        if ((m_bitmap & (1ULL << idx)) == 0ULL) {
            m_group.insert(m_group.begin() + howmany, m_group[howmany]);
            m_bitmap |= (1ULL << idx);
        }

        const int totones = bit_tools::popcnt(m_bitmap);

        // totones + 1 = dataの始点
        // m_group[howmany + 1] = idxの末尾
        m_group.insert(m_group.begin() + (totones + 1 + m_group[howmany + 1]), datvec.begin(), datvec.end());
        for (int i = howmany + 1; i <= totones; i++) {
            m_group[i] += datvec.size();
        }
    }

    std::vector<uint32_t> extract(uint32_t idx) {
        ABORT_IF_LE(SIZE, idx);

        if ((m_bitmap & (1ULL << idx)) == 0ULL) {
            return {};
        }

        const uint64_t bitmask = (1ULL << idx) - 1ULL;
        const int howmany = bit_tools::popcnt(m_bitmap & bitmask);
        const int totones = bit_tools::popcnt(m_bitmap);

        const uint32_t size = m_group[howmany + 1] - m_group[howmany];
        const uint32_t pos = totones + 1 + m_group[howmany];

        std::vector<uint32_t> vec(size);
        std::copy(&m_group[pos], &m_group[pos + size], vec.data());

        for (int i = howmany + 2; i <= totones; i++) {
            m_group[i] = m_group[i] - size;
        }
        m_group.erase(m_group.begin() + pos, m_group.begin() + pos + size);
        m_group.erase(m_group.begin() + howmany + 1);

        m_bitmap = m_bitmap & ~(1ULL << idx);

        return vec;
    }

    uint32_t size(uint32_t idx) const {
        ABORT_IF_LE(SIZE, idx);

        if ((m_bitmap & (1ULL << idx)) == 0ULL) {
            return 0;
        }

        const uint64_t bitmask = (1ULL << idx) - 1ULL;
        const int howmany = bit_tools::popcnt(m_bitmap & bitmask);

        return m_group[howmany + 1] - m_group[howmany];
    }
};

}  // namespace dyft
