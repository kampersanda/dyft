#pragma once

#include <memory>
#include <type_traits>

#include <boost/algorithm/string/classification.hpp>  // is_any_of
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>

#include "tinyformat/tinyformat.h"

#include "abort_if.hpp"
#include "splitmix.hpp"
#include "vcode_array.hpp"

namespace dyft {

inline uint64_t get_filesize(const std::string& filepath) {
    boost::filesystem::path path(filepath);
    return boost::filesystem::file_size(path);
}

inline std::ifstream make_ifstream(const std::string& filepath) {
    std::ifstream ifs(filepath);
    ABORT_IF(!ifs);
    return ifs;
}
inline std::ofstream make_ofstream(const std::string& filepath) {
    std::ofstream ofs(filepath);
    ABORT_IF(!ofs);
    return ofs;
}

inline bool exists_path(const std::string path) {
    boost::filesystem::path p(path);
    return boost::filesystem::exists(p);
}

inline void make_directory(const std::string dir) {
    boost::filesystem::path path(dir);
    if (!boost::filesystem::exists(path)) {
        if (!boost::filesystem::create_directories(path)) {
            tfm::errorfln("unable to create output directory %s", dir);
            exit(1);
        }
    }
}

template <int N>
std::unique_ptr<vcode_array<N>> load_vcodes_from_bin(const std::string& path) {
    using vint_type = typename vcode_tools<N>::vint_type;

    const uint64_t num_codes = get_filesize(path) / sizeof(uint64_t);
    std::vector<vint_type> vcodes(num_codes);

    auto ifs = make_ifstream(path);
    for (uint64_t i = 0; i < num_codes; i++) {
        uint64_t vc = 0;
        ifs.read(reinterpret_cast<char*>(&vc), sizeof(uint64_t));
        vcodes[i] = static_cast<vint_type>(vc);
    }
    return std::make_unique<vcode_array<N>>(std::move(vcodes), 1);
}

template <int N>
std::unique_ptr<vcode_array<N>> load_vcodes_from_bvecs(const std::string& path, int bits) {
    auto database = std::make_unique<vcode_array<N>>(bits);

    std::vector<uint8_t> code;
    for (auto ifs = make_ifstream(path);;) {
        uint32_t dim = 0;
        ifs.read(reinterpret_cast<char*>(&dim), sizeof(uint32_t));
        if (ifs.eof()) {
            break;
        }
        ABORT_IF_LT(dim, N);

        code.resize(dim);
        ifs.read(reinterpret_cast<char*>(code.data()), sizeof(uint8_t) * dim);

        database->append(code.data());
    }
    return database;
}

}  // namespace dyft
