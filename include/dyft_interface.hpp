#pragma once

#include "vcode_array.hpp"

namespace dyft {

template <int N>
class dyft_interface {
  public:
    using vint_type = typename vcode_tools<N>::vint_type;

    virtual ~dyft_interface() {}

    virtual uint32_t append() = 0;
    virtual void range_search(const vint_type* vcode, const std::function<void(uint32_t)>& fn) = 0;

    virtual uint32_t get_size() const = 0;
    virtual uint32_t get_leaves() const = 0;
    virtual int get_bits() const = 0;

    virtual bool selects_ls() const = 0;
    virtual double get_ls_cost() const = 0;
    virtual double get_trie_cost() const = 0;

    virtual size_t get_split_count() = 0;

    virtual void innode_stats() = 0;
    virtual void population_stats() = 0;
};

}  // namespace dyft
