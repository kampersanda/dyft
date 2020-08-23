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

    const auto codes = load_vcodes_from_bin<64>(input_path);
    std::vector<uint64_t> qcodes = sample_codes(codes->get_vcodes(), num, seed);

    auto ofs = make_ofstream(output_path);
    ofs.write(reinterpret_cast<const char*>(qcodes.data()), sizeof(uint64_t) * num);

    return 0;
}