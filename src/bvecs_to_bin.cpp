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
    p.add("input_path", "input file path of codes (in bvecs)");
    p.add("output_path", "output file path of sampled codes (in 64-binary)");
    ABORT_IF(!p.parse());

    auto input_path = p.get<std::string>("input_path");
    auto output_path = p.get<std::string>("output_path");

    const auto int_codes = load_vcodes_from_bvecs<64>(input_path, 1);

    auto ofs = make_ofstream(output_path);
    for (uint32_t i = 0; i < int_codes->get_size(); i++) {
        const uint64_t* code = int_codes->access(i);
        ofs.write(reinterpret_cast<const char*>(code), sizeof(uint64_t));
    }

    return 0;
}