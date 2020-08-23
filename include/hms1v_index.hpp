#pragma once

#include <functional>
#include <memory>
#include <unordered_map>

#include "vcode_array.hpp"

namespace dyft {

// A straightforward implementation of dynamic HmSearch 1-var (HSV), described in the paper
//  - Zhang et al. HmSearch: An efficient Hamming distance query processing algorithm, SSDBM2013
template <int N>
class hms1v_index {
  public:
    struct idvec_pair {
        std::vector<uint32_t> exact_vec;
        std::vector<uint32_t> variant_vec;
    };
    using vint_type = typename vcode_tools<N>::vint_type;
    using table_type = std::unordered_map<uint64_t, idvec_pair>;

  private:
    const vcode_array<N>* m_database = nullptr;
    const int m_sigma = 0;  // deletion marker also
    const int m_radius = 0;
    const int m_blocks = 0;

    std::vector<table_type> m_tables;
    std::vector<int> m_begs;
    uint32_t m_ids = 0;

    std::unordered_map<uint32_t, std::vector<uint32_t>> m_cand_map;

  public:
    hms1v_index(const vcode_array<N>* database, int radius)
        : m_database(database), m_sigma(1 << m_database->get_bits()), m_radius(radius), m_blocks((radius + 3) / 2),
          m_tables(m_blocks), m_begs(m_blocks + 1) {
        m_begs[0] = 0;
        for (int b = 0; b < m_blocks; b++) {
            const int m = (b + N) / m_blocks;
            m_begs[b + 1] = m_begs[b] + m;
        }
        ABORT_IF_NE(m_begs[m_blocks], N);
    }

    uint32_t append() {
        const uint32_t new_id = m_ids++;
        const vint_type* vcode = m_database->access(new_id);

        int code[N];
        for (int i = 0; i < N; i++) {
            code[i] = vcode_tools<N>::get_int(vcode, i, get_bits());
        }
        for (int b = 0; b < m_blocks; b++) {
            insert(code, b, new_id);
        }

        return new_id;
    }

    void range_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) {
        int code[N];
        for (int i = 0; i < N; i++) {
            code[i] = vcode_tools<N>::get_int(vcode, i, get_bits());
        }

        m_cand_map.clear();

        for (int b = 0; b < m_blocks; b++) {
            search(code, b, [&](uint32_t id, uint32_t errs) {
                auto it = m_cand_map.find(id);
                if (it != m_cand_map.end()) {
                    it->second.push_back(errs);
                } else {
                    m_cand_map.insert(std::make_pair(id, std::vector<uint32_t>{errs}));
                }
            });
        }

        for (const auto& kv : m_cand_map) {
            uint32_t cand_id = kv.first;
            const std::vector<uint32_t>& errors = kv.second;

            ABORT_IF_EQ(errors.size(), 0);

            // enhanced filter
            bool filtered = false;

            if (m_radius % 2 == 0) {
                if (errors.size() < 2) {  // has less than two number
                    if (errors[0] == 1) {
                        filtered = true;
                    }
                }
            } else {
                if (errors.size() < 3) {  // has less than three number
                    if (errors.size() == 1) {
                        filtered = true;
                    } else if (errors[0] == 1 and errors[1] == 1) {
                        filtered = true;
                    }
                }
            }

            if (!filtered) {
                fn(cand_id);
            }
        }
    }

    uint32_t get_size() const {
        return m_ids;
    }
    int get_bits() const {
        return m_database->get_bits();
    }

  private:
    static uint64_t fnv1a(const int* str, const int n) {
        static const uint64_t init = uint64_t((sizeof(uint64_t) == 8) ? 0xcbf29ce484222325ULL : 0x811c9dc5ULL);
        static const uint64_t multiplier = uint64_t((sizeof(uint64_t) == 8) ? 0x100000001b3ULL : 0x1000193ULL);

        uint64_t h = init;
        for (int i = 0; i < n; i++) {
            h ^= uint64_t(str[i]);
            h *= multiplier;
        }
        return h;
    }

    void insert(const int* code, const int bpos, const uint32_t new_id) {
        const int beg = m_begs[bpos];
        const int len = m_begs[bpos + 1] - m_begs[bpos];

        code += beg;
        auto& table = m_tables[bpos];

        // exact one
        {
            const uint64_t h = fnv1a(code, len);

            auto it = table.find(h);
            if (it != table.end()) {  // found?
                it->second.exact_vec.push_back(new_id);
            } else {
                table.insert(std::make_pair(h, idvec_pair{std::vector<uint32_t>{new_id}, std::vector<uint32_t>{}}));
            }
        }

        // 1-variant ones
        int sig[N];
        for (int i = 0; i < len; i++) {
            std::copy(code, code + len, sig);

            for (int c = 0; c < m_sigma; c++) {
                if (c == code[i]) {
                    continue;
                }

                sig[i] = c;
                const uint64_t h = fnv1a(sig, len);

                auto it = table.find(h);
                if (it != table.end()) {  // found?
                    it->second.variant_vec.push_back(new_id);
                } else {
                    table.insert(std::make_pair(h, idvec_pair{std::vector<uint32_t>{}, std::vector<uint32_t>{new_id}}));
                }
            }
        }
    }

    void search(const int* code, const int bpos, const std::function<void(uint32_t, bool)>& fn) {
        const int beg = m_begs[bpos];
        const int len = m_begs[bpos + 1] - m_begs[bpos];

        code += beg;
        const auto& table = m_tables[bpos];

        const uint64_t h = fnv1a(code, len);

        auto it = table.find(h);
        if (it == table.end()) {  // not found?
            return;
        }

        for (uint32_t id : it->second.exact_vec) {
            fn(id, 0);
        }
        for (uint32_t id : it->second.variant_vec) {
            fn(id, 1);
        }
    }
};

}  // namespace dyft
