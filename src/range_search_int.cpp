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

static constexpr double ABORT_THRESHOLD_MS = 100.0;

constexpr uint32_t SCALES[] = {
    10, 100, 1'000, 10'000, 100'000, 1'000'000, 10'000'000, 100'000'000, 1'000'000'000, UINT32_MAX,
};

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("base_path", "input file path of database (in 64-bit binary)");
    p.add("query_path", "input file path of queries (in 64-bit binary)");
    p.add("result_dir", "output directory path of results", "-o", false);
    p.add("radius", "Hamming distance threshold", "-R", false);
    p.add("ints", "number of ints for each vector (32|64)", "-N", false);
    p.add("bits", "number of bits for each integer [1,8]", "-B", false);
#if defined(ALGO_DYFT)
    p.add("algorithm", "algorithm (array|art|mart)", "-A", false);
    p.add("splitthr", "split threshold (0 means optimal assignment)", "-T", false);
    p.add("in_weight", "weigth of innode", "-W", false);
    p.add("blocks", "number of blocks", "-K", false);
#endif
    ABORT_IF(!p.parse());
    return p;
}

template <int N>
int main_template(const cmd_line_parser::parser& p) {
    using vint_type = typename vcode_traits<N>::vint_type;

    const auto base_path = p.get<std::string>("base_path");
    const auto query_path = p.get<std::string>("query_path");
    const auto result_dir = p.get<std::string>("result_dir", "results");
    const auto radius = p.get<int>("radius", 2);
    const auto bits = p.get<int>("bits", 4);
#if defined(ALGO_DYFT)
    const auto algorithm = p.get<std::string>("algorithm", "mart");
    const auto splitthr = p.get<int>("splitthr", 0);
    const auto in_weight = p.get<double>("in_weight", 1.0);
    const auto blocks = p.get<int>("blocks", 2);
#endif

#ifdef ALGO_LS
    tfm::printfln("algorithm: ls");
    STATISTIC_TAG("algorithm", "ls");
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
    const auto query_codes = load_vcodes_from_bvecs<N>(query_path, bits);
    tfm::printfln("base_path: %s", base_path);
    tfm::printfln("query_path: %s", query_path);
    STATISTIC_TAG("base_path", base_path);
    STATISTIC_TAG("query_path", query_path);

    const uint32_t base_size = base_codes->get_size();
    const uint32_t query_size = query_codes->get_size();
    tfm::printfln("base_size: %d", base_size);
    tfm::printfln("query_size: %d", query_size);
    STATISTIC_TAG("base_size", base_size);
    STATISTIC_TAG("query_size", query_size);

    std::vector<uint32_t> test_scales;
    for (uint32_t i = 0; SCALES[i] < base_size; i++) {
        test_scales.push_back(SCALES[i]);
    }
    test_scales.push_back(base_size);

    uint32_t prev_size = 0;
    std::vector<uint32_t> counts(query_size);
    std::vector<uint32_t> verify_counts(query_size);

#ifdef ALGO_LS
    double search_time_in_ms = 0.0;

    for (uint32_t test_size : test_scales) {
        tfm::printfln("# %d codes...", test_size);

        START_TIMER(Search);
        for (uint32_t qi = 0; qi < query_size; qi++) {
            const vint_type* qvcode = query_codes->access(qi);
            for (uint32_t bi = prev_size; bi < test_size; bi++) {
                const vint_type* bvcode = base_codes->access(bi);
                const int hamdist = vcode_tools<N>::get_hamdist(qvcode, bvcode, bits, radius);
                if (hamdist <= radius) {
                    counts[qi] += 1;
                }
            }
        }
        STOP_TIMER_V(Search);

        search_time_in_ms += GET_TIMER_MILLISEC(Search);
        const double ave_search_time_in_ms = search_time_in_ms / query_size;

        const double ave_num_results = get_average(counts);
        tfm::reportfln("average number of results: %g", ave_num_results);
        tfm::reportfln("average search time in ms: %g", ave_search_time_in_ms);

        STATISTIC_APPEND("search", {{"num_codes", test_size},  //
                                    {"ave_search_time_in_ms", ave_search_time_in_ms},  //
                                    {"ave_num_results", ave_num_results}});

        if (ave_search_time_in_ms >= ABORT_THRESHOLD_MS) {
            tfm::warnfln("Abort search becasue the time exceeds %g ms", ABORT_THRESHOLD_MS);
            break;
        }

        prev_size = test_size;
    }

    const auto result_path = tfm::format("%s/range_search_int_ls-%s-%s-%dN-%dB-%dR.json",  //
                                         result_dir, normalize_filepath(base_path), normalize_filepath(query_path),  //
                                         N, bits, radius);
#endif

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

    double insert_time_in_ms = 0.0;

    for (uint32_t test_size : test_scales) {
        tfm::printfln("# %d codes...", test_size);

        START_TIMER(Insert);
        for (uint32_t bi = prev_size; bi < test_size; bi++) {
            ABORT_IF_NE(bi, index->append());
        }
        STOP_TIMER_V(Insert);

        insert_time_in_ms += GET_TIMER_MILLISEC(Insert);
        const double ave_insert_time_in_ms = insert_time_in_ms / test_size;

        START_TIMER(Search);
        for (uint32_t qi = 0; qi < query_size; qi++) {
            const vint_type* qvcode = query_codes->access(qi);

            uint32_t count = 0;
            uint32_t verify_count = 0;

            index->range_search(qvcode, [&](uint32_t bi) {
                const vint_type* bvcode = base_codes->access(bi);
                const int hamdist = vcode_tools<N>::get_hamdist(qvcode, bvcode, bits, radius);
                if (hamdist <= radius) {
                    count += 1;
                }
                verify_count += 1;
            });

            counts[qi] = count;
            verify_counts[qi] = verify_count;
        }
        STOP_TIMER_V(Search);

        std::string selected = index->selects_ls() ? "LS" : "Trie";
        tfm::reportfln("selected search: %s", selected);

        tfm::reportfln("ls_cost:   %g", index->get_ls_cost());
        tfm::reportfln("trie_cost: %g", index->get_trie_cost());

        const double ave_search_time_in_ms = GET_TIMER_MILLISEC(Search) / query_size;
        tfm::reportfln("average insert time in ms: %g", ave_insert_time_in_ms);
        tfm::reportfln("average search time in ms: %g", ave_search_time_in_ms);

        const double ave_num_results = get_average(counts);
        const double ave_num_cadidates = get_average(verify_counts);
        tfm::reportfln("average number of results: %g", ave_num_results);
        tfm::reportfln("average number of candidates: %g", ave_num_cadidates);

        tfm::reportfln("num leaves: %d", index->get_leaves());
        tfm::reportfln("num splits: %d", index->get_split_count());

        STATISTIC_APPEND("search", {{"num_codes", test_size},  //
                                    {"num_leaves", index->get_leaves()},  //
                                    {"num_splits", index->get_split_count()},  //
                                    {"selected_search", selected},  //
                                    {"ave_insert_time_in_ms", ave_insert_time_in_ms},  //
                                    {"ave_search_time_in_ms", ave_search_time_in_ms},  //
                                    {"ave_num_results", ave_num_results},  //
                                    {"ave_num_cadidates", ave_num_cadidates}});

        if (ave_search_time_in_ms >= ABORT_THRESHOLD_MS) {
            tfm::warnfln("Abort search becasue the time exceeds %g ms", ABORT_THRESHOLD_MS);
            break;
        }

        prev_size = test_size;
    }

    const auto result_path = tfm::format("%s/range_search_int_%s-%s-%s-%dN-%dB-%dR-%dT-%gW-%dK.json",  //
                                         result_dir, algorithm,  //
                                         normalize_filepath(base_path), normalize_filepath(query_path),  //
                                         N, bits, radius, splitthr, in_weight * 100, blocks);
#endif  // defined(ALGO_DYFT)

#if defined(ALGO_HMS1DV)
    try {
        double insert_time_in_ms = 0.0;
        hms1dv_index<N> index(base_codes.get(), radius);

        for (uint32_t test_size : test_scales) {
            tfm::printfln("# %d codes...", test_size);

            START_TIMER(Insert);
            for (uint32_t bi = prev_size; bi < test_size; bi++) {
                ABORT_IF_NE(bi, index.append());
            }
            STOP_TIMER_V(Insert);

            insert_time_in_ms += GET_TIMER_MILLISEC(Insert);
            const double ave_insert_time_in_ms = insert_time_in_ms / test_size;

            START_TIMER(Search);
            for (uint32_t qi = 0; qi < query_size; qi++) {
                const vint_type* qvcode = query_codes->access(qi);

                uint32_t count = 0;
                uint32_t verify_count = 0;

                index.range_search(qvcode, [&](uint32_t bi) {
                    const vint_type* bvcode = base_codes->access(bi);
                    const int hamdist = vcode_tools<N>::get_hamdist(qvcode, bvcode, bits, radius);
                    if (hamdist <= radius) {
                        count += 1;
                    }
                    verify_count += 1;
                });
                counts[qi] = count;
                verify_counts[qi] = verify_count;
            }
            STOP_TIMER_V(Search);

            const double ave_search_time_in_ms = GET_TIMER_MILLISEC(Search) / query_size;
            tfm::reportfln("average insert time in ms: %g", ave_insert_time_in_ms);
            tfm::reportfln("average search time in ms: %g", ave_search_time_in_ms);

            const double ave_num_results = get_average(counts);
            const double ave_num_cadidates = get_average(verify_counts);
            tfm::reportfln("average number of results: %g", ave_num_results);
            tfm::reportfln("average number of candidates: %g", ave_num_cadidates);

            STATISTIC_APPEND("search", {{"num_codes", test_size},  //
                                        {"ave_insert_time_in_ms", ave_insert_time_in_ms},  //
                                        {"ave_search_time_in_ms", ave_search_time_in_ms},  //
                                        {"ave_num_results", ave_num_results},  //
                                        {"ave_num_cadidates", ave_num_cadidates}});

            if (ave_search_time_in_ms >= ABORT_THRESHOLD_MS) {
                tfm::warnfln("Abort search becasue the time exceeds %g ms", ABORT_THRESHOLD_MS);
                break;
            }

            prev_size = test_size;
        }
    } catch (const std::bad_alloc& e) {
        tfm::errorfln("Allocation failed: %s", e.what());
    }

    const auto result_path = tfm::format("%s/range_search_int_hms1dv-%s-%s-%dN-%dB-%dR.json",  //
                                         result_dir, normalize_filepath(base_path), normalize_filepath(query_path),  //
                                         N, bits, radius);
#endif  // defined(ALGO_HMS1DV)

#if defined(ALGO_GV)
    double insert_time_in_ms = 0.0;
    gv_index<N> index(base_codes.get(), radius);

    for (uint32_t test_size : test_scales) {
        tfm::printfln("# %d codes...", test_size);

        START_TIMER(Insert);
        for (uint32_t bi = prev_size; bi < test_size; bi++) {
            ABORT_IF_NE(bi, index.append());
        }
        STOP_TIMER_V(Insert);

        insert_time_in_ms += GET_TIMER_MILLISEC(Insert);
        const double ave_insert_time_in_ms = insert_time_in_ms / test_size;

        START_TIMER(Search);
        for (uint32_t qi = 0; qi < query_size; qi++) {
            const vint_type* qvcode = query_codes->access(qi);

            uint32_t count = 0;
            uint32_t verify_count = 0;

            index.range_search(qvcode, [&](uint32_t bi) {
                const vint_type* bvcode = base_codes->access(bi);
                const int hamdist = vcode_tools<N>::get_hamdist(qvcode, bvcode, bits, radius);
                if (hamdist <= radius) {
                    count += 1;
                }
                verify_count += 1;
            });
            counts[qi] = count;
            verify_counts[qi] = verify_count;
        }
        STOP_TIMER_V(Search);

        const double ave_search_time_in_ms = GET_TIMER_MILLISEC(Search) / query_size;
        tfm::reportfln("average insert time in ms: %g", ave_insert_time_in_ms);
        tfm::reportfln("average search time in ms: %g", ave_search_time_in_ms);

        const double ave_num_results = get_average(counts);
        const double ave_num_cadidates = get_average(verify_counts);
        tfm::reportfln("average number of results: %g", ave_num_results);
        tfm::reportfln("average number of candidates: %g", ave_num_cadidates);

        STATISTIC_APPEND("search", {{"num_codes", test_size},  //
                                    {"ave_insert_time_in_ms", ave_insert_time_in_ms},  //
                                    {"ave_search_time_in_ms", ave_search_time_in_ms},  //
                                    {"ave_num_results", ave_num_results},  //
                                    {"ave_num_cadidates", ave_num_cadidates}});

        if (ave_search_time_in_ms >= ABORT_THRESHOLD_MS) {
            tfm::warnfln("Abort search becasue the time exceeds %g ms", ABORT_THRESHOLD_MS);
            break;
        }

        prev_size = test_size;
    }

    const auto result_path = tfm::format("%s/range_search_int_gv-%s-%s-%dN-%dB-%dR.json",  //
                                         result_dir, normalize_filepath(base_path), normalize_filepath(query_path),  //
                                         N, bits, radius);
#endif  // defined(ALGO_GV)

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