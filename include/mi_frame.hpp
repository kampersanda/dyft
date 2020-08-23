#pragma once

#include "array_index.hpp"
#include "art_index.hpp"
#include "mart_index.hpp"

namespace dyft {

template <class Index>
class mi_frame : public dyft_interface<Index::LEN> {
  public:
    using index_type = Index;
    using vint_type = typename index_type::vint_type;

    static constexpr int LEN = index_type::LEN;

  private:
    const vcode_array<LEN>* m_database = nullptr;
    const int m_blocks = 0;

    std::vector<int> m_radii;
    std::vector<std::unique_ptr<index_type>> m_indexes;
    uint32_t m_ids = 0;

    // for query
    std::vector<uint32_t> m_cands;

  public:
    mi_frame(const vcode_array<LEN>* database, int radius, int blocks, int splitthr, double in_weight)
        : m_database(database), m_blocks(blocks), m_radii(blocks), m_indexes(blocks) {
        ABORT_IF_LE(blocks, 1);

        int beg = 0;  // of chunk positios
        const int gph = radius - m_blocks + 1;
        const int chunks_size = index_type::get_chunks_size(m_database->get_bits());

        for (int b = 0; b < m_blocks; b++) {
            const int m = (b + chunks_size) / m_blocks;
            m_radii[b] = (b + gph) / m_blocks;
            if (m_radii[b] >= 0) {
                m_indexes[b] = std::make_unique<index_type>(database, m_radii[b], beg, beg + m, splitthr, in_weight);
            }
            beg += m;
        }
        ABORT_IF_NE(beg, chunks_size);

        m_cands.reserve(1U << 10);
    }

    uint32_t append() override {
        for (int b = 0; b < m_blocks; b++) {
            if (m_indexes[b]) {
                ABORT_IF_NE(m_indexes[b]->append(), m_ids);
            }
        }
        return m_ids++;
    }

    void range_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) override {
        m_cands.clear();

        if (selects_ls()) {
            for (uint32_t id = 0; id < get_size(); id++) {
                fn(id);
            }
            return;
        }

        for (int b = 0; b < m_blocks; b++) {
            if (m_radii[b] >= 0) {
                DEBUG_ABORT_IF(!m_indexes[b]);
                m_indexes[b]->trie_search(vcode, [&](uint32_t bi) { m_cands.push_back(bi); });
            }
        }

        std::sort(m_cands.begin(), m_cands.end());
        m_cands.erase(std::unique(m_cands.begin(), m_cands.end()), m_cands.end());

        for (uint32_t id : m_cands) {
            fn(id);
        }
    }

    uint32_t get_size() const override {
        return m_ids;
    }
    int get_bits() const override {
        return m_database->get_bits();
    }

    uint32_t get_leaves() const override {
        return 0;
    }
    size_t get_split_count() override {
        return 0;
    }

    bool selects_ls() const override {
        for (int b = 0; b < m_blocks; b++) {
            if (m_indexes[b] and m_indexes[b]->selects_ls()) {
                return true;
            }
        }
        return false;
    }

    double get_ls_cost() const override {
        return 0.0;
    }
    double get_trie_cost() const override {
        return 0.0;
    }

    void innode_stats() override {}
    void population_stats() override {}
};

}  // namespace dyft
