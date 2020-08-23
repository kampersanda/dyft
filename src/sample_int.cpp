#include <random>

#include <cmd_line_parser/parser.hpp>

#include <abort_if.hpp>
#include <io.hpp>
#include <vcode_array.hpp>

using namespace dyft;

std::vector<uint64_t> sample_codes(const std::vector<uint64_t>& codes, uint32_t num, uint64_t seed) {
    splitmix64 engine(seed);
    std::vector<uint64_t> sampled(num);
    for (uint32_t i = 0; i < num; ++i) {
        sampled[i] = codes[engine.next() % codes.size()];
    }
    return sampled;
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "input file path of simhashed codes");
    p.add("output_path", "output file path of sampled simhashed codes");
    p.add("num", "number of sample codes");
    p.add("seed", "seed", "-s", false);
    ABORT_IF(!p.parse());

    auto input_path = p.get<std::string>("input_path");
    auto output_path = p.get<std::string>("output_path");
    auto num = p.get<size_t>("num");
    auto seed = p.get<size_t>("seed", 114514);

    std::vector<uint8_t> codes;
    uint32_t size = 0;
    uint32_t dim = 0;

    for (auto ifs = make_ifstream(input_path);;) {
        ifs.read(reinterpret_cast<char*>(&dim), sizeof(uint32_t));
        if (ifs.eof()) {
            break;
        }
        ABORT_IF_LT(256, dim);

        uint8_t code[256];
        ifs.read(reinterpret_cast<char*>(code), sizeof(uint8_t) * dim);
        std::copy(code, code + dim, std::back_inserter(codes));

        size += 1;
    }

    std::default_random_engine engine(seed);
    std::uniform_int_distribution<size_t> dist(0, size - 1);

    auto ofs = make_ofstream(output_path);
    for (size_t i = 0; i < num; ++i) {
        const uint64_t pos = dist(engine) * dim;
        ofs.write(reinterpret_cast<const char*>(&dim), sizeof(uint32_t));
        ofs.write(reinterpret_cast<const char*>(codes.data() + pos), sizeof(uint8_t) * dim);
    }

    return 0;
}