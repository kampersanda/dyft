#include <tinyformat/tinyformat.h>
#include <cmd_line_parser/parser.hpp>

#include <io.hpp>
#include <misc.hpp>
#include <statistic_reporter.hpp>
#include <timer.hpp>

#ifdef ALGO_DYFT
#include <mi_frame.hpp>
#endif
#ifdef ALGO_HMS1DV
#include <hms1dv_index.hpp>
#endif
#ifdef ALGO_GV
#include <gv_index.hpp>
#endif

using namespace dyft;

static constexpr size_t ABORT_THRESHOLD_GiB = 100;

constexpr uint32_t SCALES[] = {
    10, 100, 1'000, 10'000, 100'000, 1'000'000, 10'000'000, 100'000'000, 1'000'000'000, UINT32_MAX,
};

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("base_path", "input file path of database (in 64-bit binary)");
    p.add("result_dir", "output directory path of results", "-o", false);
    p.add("radius", "Hamming distance threshold", "-R", false);
    p.add("ints", "number of ints for each vector (32|64)", "-N", false);
    p.add("bits", "number of bits for each integer [1,8]", "-B", false);
#if defined(ALGO_DYFT)
    p.add("algorithm", "algorithm (array|art|mart)", "-A", false);
    p.add("splitthr", "split threshold (0 means optimal assignment)", "-T", false);
    p.add("in_weight", "weigth of innode", "-W", false);
    p.add("blocks", "number of blocks (0 means reasonable number based on radius)", "-K", false);
#endif
    ABORT_IF(!p.parse());
    return p;
}

template <int N>
int main_template(const cmd_line_parser::parser& p) {
    const auto base_path = p.get<std::string>("base_path");
    const auto result_dir = p.get<std::string>("result_dir", "results");
    const auto radius = p.get<int>("radius", 2);
    const auto bits = p.get<int>("bits", 4);
#if defined(ALGO_DYFT)
    const auto algorithm = p.get<std::string>("algorithm", "mart");
    const auto splitthr = p.get<int>("splitthr", 0);
    const auto in_weight = p.get<double>("in_weight", 1.0);
    const auto blocks = p.get<int>("blocks", 2);
#endif

#ifdef ALGO_DYFT
    ABORT_IF((algorithm != "array") and (algorithm != "art") and (algorithm != "mart"));
    tfm::printfln("algorithm: %s", algorithm);
    tfm::printfln("splitthr: %d", splitthr);
    tfm::printfln("in_weight: %g", in_weight);
    tfm::printfln("blocks: %d", blocks);
    STATISTIC_TAG("algorithm", algorithm);
    STATISTIC_TAG("splitthr", splitthr);
    STATISTIC_TAG("in_weight", in_weight);
    STATISTIC_TAG("blocks", blocks);
#endif
#ifdef ALGO_HMS1DV
    tfm::printfln("algorithm: hms1dv");
    STATISTIC_TAG("algorithm", "hms1dv");
#endif
#ifdef ALGO_GV
    tfm::printfln("algorithm: gv");
    STATISTIC_TAG("algorithm", "gv");
#endif

    tfm::printfln("ints: %d", N);
    tfm::printfln("bits: %d", bits);
    tfm::printfln("radius: %d", radius);
    STATISTIC_TAG("ints", N);
    STATISTIC_TAG("bits", bits);
    STATISTIC_TAG("radius", radius);

    const auto base_codes = load_vcodes_from_bvecs<N>(base_path, bits);
    tfm::printfln("base_path: %s", base_path);
    STATISTIC_TAG("base_path", base_path);

    const uint32_t base_size = base_codes->get_size();
    tfm::printfln("base_size: %d", base_size);
    STATISTIC_TAG("base_size", base_size);

    std::vector<uint32_t> test_scales;
    for (uint32_t i = 0; SCALES[i] < base_size; i++) {
        test_scales.push_back(SCALES[i]);
    }
    test_scales.push_back(base_size);

    std::vector<uint64_t> process_sizes(test_scales.size());
    std::vector<double> insertion_times(test_scales.size());

    uint32_t prev_size = 0;
    uint64_t init_process_bytes = get_process_size_in_bytes();
    double insert_time_in_sec = 0.0;

#ifdef ALGO_DYFT
    std::unique_ptr<dyft_interface<N>> index;

    if (algorithm == "array") {
        if (blocks == 1) {
            index = std::make_unique<array_index<N>>(base_codes.get(), radius, splitthr, in_weight);
        } else {
            index = std::make_unique<mi_frame<array_index<N>>>(
                base_codes.get(), radius, blocks == 0 ? radius / 2 + 1 : blocks, splitthr, in_weight);
        }
    } else if (algorithm == "art") {
        if (blocks == 1) {
            index = std::make_unique<art_index<N>>(base_codes.get(), radius, splitthr, in_weight);
        } else {
            index = std::make_unique<mi_frame<art_index<N>>>(
                base_codes.get(), radius, blocks == 0 ? radius / 2 + 1 : blocks, splitthr, in_weight);
        }
    } else {  // algorithm == "mart"
        if (blocks == 1) {
            index = std::make_unique<mart_index<N>>(base_codes.get(), radius, splitthr, in_weight);
        } else {
            index = std::make_unique<mi_frame<mart_index<N>>>(
                base_codes.get(), radius, blocks == 0 ? radius / 2 + 1 : blocks, splitthr, in_weight);
        }
    }
#endif
#if defined(ALGO_HMS1DV)
    auto index = std::make_unique<hms1dv_index<N>>(base_codes.get(), radius);
#endif
#if defined(ALGO_GV)
    auto index = std::make_unique<gv_index<N>>(base_codes.get(), radius);
#endif

    for (uint32_t i = 0; i < test_scales.size(); i++) {
        const uint32_t test_size = test_scales[i];
        tfm::printfln("# %d codes...", test_size);

        START_TIMER(Insert);
        for (uint32_t bi = prev_size; bi < test_size; bi++) {
            ABORT_IF_NE(bi, index->append());
        }
        STOP_TIMER_V(Insert);

        const size_t process_size = get_process_size_in_bytes();
        process_sizes[i] = process_size - init_process_bytes;

        insert_time_in_sec += GET_TIMER_SEC(Insert);
        insertion_times[i] = insert_time_in_sec;

        tfm::reportfln("process size in MiB: %g", process_sizes[i] / (1024.0 * 1024.0));
        tfm::reportfln("insertion time in sec: %g", insertion_times[i]);

        if (to_GiB(process_size) >= ABORT_THRESHOLD_GiB) {
            tfm::warnfln("Abort build becasue the memory exceeds %d GiB", ABORT_THRESHOLD_GiB);
            test_scales.resize(i + 1);
            break;
        }

        prev_size = test_size;
    }

    for (uint32_t i = 0; i < test_scales.size(); i++) {
        const uint32_t test_size = test_scales[i];
        const uint64_t process_size = process_sizes[i];
        const double insertion_time = insertion_times[i];
        STATISTIC_APPEND("insertion", {{"num_codes", test_size},  //
                                       {"process_bytes", process_size},  //
                                       {"insertion_time_in_sec", insertion_time}});
    }

#ifdef ALGO_DYFT
    const auto result_path = tfm::format("%s/build_index_int_%s-%s-%dN-%dB-%dR-%dT-%gW-%dK.json",  //
                                         result_dir, algorithm, normalize_filepath(base_path),  //
                                         N, bits, radius, splitthr, in_weight * 100, blocks);
#endif
#ifdef ALGO_HMS1DV
    const auto result_path = tfm::format("%s/build_index_int_hms1dv-%s-%dN-%dB-%dR.json",  //
                                         result_dir, normalize_filepath(base_path), N, bits, radius);
#endif
#ifdef ALGO_GV
    const auto result_path = tfm::format("%s/build_index_int_gv-%s-%dN-%dB-%dR.json",  //
                                         result_dir, normalize_filepath(base_path), N, bits, radius);
#endif

    make_directory(result_dir);
    STATISTIC_SAVE(result_path);
    tfm::printfln("wrote %s", result_path);

    return 0;
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    auto p = make_parser(argc, argv);
    const auto ints = p.get<int>("ints", 64);

    switch (ints) {
        case 32:
            return main_template<32>(p);
        case 64:
            return main_template<64>(p);
        default:
            break;
    }

    return 1;
}