#pragma once

#include <functional>
#include <memory>
#include <unordered_map>

#include "vcode_array.hpp"

namespace dyft {

// A straightforward implementation of multi-index hashing for integer sketches, following the idea
//  - Gog and Venturini. Fast and compact Hamming distance index, SIGIR16
template <int N>
class gv_index {
  public:
    using vint_type = typename vcode_tools<N>::vint_type;
    using table_type = std::unordered_map<uint64_t, std::vector<uint32_t>>;

  private:
    const vcode_array<N>* m_database = nullptr;
    const int m_sigma = 0;
    const int m_radius = 0;
    const int m_blocks = 0;

    std::vector<table_type> m_tables;
    std::vector<int> m_begs;
    uint32_t m_ids = 0;

  public:
    gv_index(const vcode_array<N>* database, int radius)
        : m_database(database), m_sigma(1 << m_database->get_bits()), m_radius(radius), m_blocks(radius / 2 + 1),
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
        std::vector<uint32_t> cands;
        cands.reserve(1U << 10);

        int code[N];
        for (int i = 0; i < N; i++) {
            code[i] = vcode_tools<N>::get_int(vcode, i, get_bits());
        }
        for (int b = 0; b < m_blocks; b++) {
            search(code, b, [&](uint32_t id) { cands.push_back(id); });
        }

        std::sort(cands.begin(), cands.end());
        cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

        for (uint32_t id : cands) {
            fn(id);
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

        const uint64_t h = fnv1a(code, len);

        auto it = table.find(h);
        if (it != table.end()) {  // found?
            it->second.push_back(new_id);
        } else {
            table.insert(std::make_pair(h, std::vector<uint32_t>{new_id}));
        }
    }

    void search(const int* code, const int bpos, const std::function<void(uint32_t)>& fn) {
        const int beg = m_begs[bpos];
        const int len = m_begs[bpos + 1] - m_begs[bpos];

        code += beg;
        const auto& table = m_tables[bpos];

        // Exact Search
        {
            const uint64_t h = fnv1a(code, len);
            auto it = table.find(h);
            if (it != table.end()) {  // found?
                for (uint32_t id : it->second) {
                    fn(id);
                }
            }
        }

        // 1-err searches
        int sig[N];
        for (int i = 0; i < len; i++) {
            std::copy(code, code + len, sig);
            for (int j = 1; j < m_sigma; j++) {
                sig[i] = (sig[i] + 1) % m_sigma;
                const uint64_t h = fnv1a(sig, len);
                auto it = table.find(h);
                if (it != table.end()) {  // found?
                    for (uint32_t id : it->second) {
                        fn(id);
                    }
                }
            }
        }
    }
};

}  // namespace dyft
