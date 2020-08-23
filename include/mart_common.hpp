#pragma once

#include "abort_if.hpp"

namespace dyft {

/* * * * * * * * * * * *
 *  Basic definitions of MART
 */
enum mart_insert_flags {
    MART_FOUND,
    MART_INSERTED,
    MART_NEEDED_TO_EXPAND,
};

#define MART_SPACE_EFFICIENT

#ifdef MART_SPACE_EFFICIENT
enum mart_node_types : uint8_t {
    MART_LEAF,
    MART_NODE_2,
    MART_NODE_4,
    MART_NODE_8,
    MART_NODE_16,
    MART_NODE_32,
    MART_NODE_64,
    MART_NODE_128,
    MART_NODE_256,
    MART_NIL_NTYPE,
};
#else
enum mart_node_types : uint8_t {
    MART_LEAF,
    MART_NODE_4,
    MART_NODE_32,
    MART_NODE_256,
    MART_NIL_NTYPE,
};
#endif

static constexpr uint32_t MART_NILID = UINT32_MAX;

struct mart_pointer {
    uint32_t nid;
    mart_node_types ntype;
};

inline mart_pointer make_mart_nullptr() {
    return mart_pointer{MART_NILID, MART_NIL_NTYPE};
}
inline bool check_mart_nullptr(mart_pointer p) {
    return p.nid == MART_NILID;
}

struct mart_cursor {
    uint32_t offset;
    mart_pointer pptr;  // src pointer
    mart_pointer nptr;
};

struct mart_edge {
    uint8_t c;
    mart_pointer nptr;
};

}  // namespace dyft
