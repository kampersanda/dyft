#pragma once

#include <array>
#include <unordered_map>

#include "abort_if.hpp"
#include "vcode_tools.hpp"

namespace dyft {

template <int N>
class vcode_array {
  public:
    using vint_type = typename vcode_tools<N>::vint_type;

  private:
    uint32_t m_size = 0;
    std::vector<vint_type> m_vcodes;

    const int m_bits = 0;

  public:
    explicit vcode_array(int bits) : m_bits(bits) {
        ABORT_IF_OUT(m_bits, 1, 8);
    }
    explicit vcode_array(std::vector<vint_type>&& vcodes, int bits)
        : m_size(vcodes.size() / bits), m_vcodes(std::move(vcodes)), m_bits(bits) {}

    uint32_t append(const uint8_t* code) {
        vint_type vcode[8];
        vcode_tools<N>::to_vints(code, vcode, m_bits);
        std::copy(vcode, vcode + m_bits, std::back_inserter(m_vcodes));
        return m_size++;
    }

    const vint_type* access(uint32_t id) const {
        ABORT_IF_LE(m_size, id);
        return &m_vcodes[id * m_bits];
    }

    const vint_type* data() const {
        return m_vcodes.data();
    }

    uint32_t get_size() const {
        return m_size;
    }
    int get_bits() const {
        return m_bits;
    }

    const std::vector<vint_type>& get_vcodes() const {
        return m_vcodes;
    }
};

}  // namespace dyft
