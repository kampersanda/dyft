#include <random>

#include <cmd_line_parser/parser.hpp>

#include <abort_if.hpp>
#include <io.hpp>
#include <vcode_array.hpp>

using namespace dyft;

static constexpr uint32_t DIM = 64;

int main(int argc, char** argv) {
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    cmd_line_parser::parser p(argc, argv);
    p.add("output_path", "output file path of random codes (in bvecs)");
    p.add("num", "number of codes");
    p.add("seed", "seed", "-s", false);
    ABORT_IF(!p.parse());

    auto output_path = p.get<std::string>("output_path");
    auto num = p.get<size_t>("num");
    auto seed = p.get<size_t>("seed", 114514);

    std::default_random_engine engine(seed);
    std::uniform_int_distribution<> dist(0, 255);

    auto ofs = make_ofstream(output_path);

    for (size_t i = 0; i < num; i++) {
        uint8_t code[DIM];
        for (size_t j = 0; j < DIM; j++) {
            code[j] = static_cast<uint8_t>(dist(engine));
        }
        ofs.write(reinterpret_cast<const char*>(&DIM), sizeof(uint32_t));
        ofs.write(reinterpret_cast<const char*>(code), sizeof(uint8_t) * DIM);
    }

    return 0;
}