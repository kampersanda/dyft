#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>

#include <tinyformat/tinyformat.h>
#include <cmd_line_parser/parser.hpp>

#include <abort_if.hpp>
#include <io.hpp>

using namespace dyft;

const size_t REPETITIONS = 64;  // i.e., output code length

class fvecs_iterator {
  public:
    fvecs_iterator(const std::string& path, size_t dimension) : m_ifs(make_ifstream(path)), m_vec(dimension) {}

    const std::vector<float>& next() {
        uint32_t dim = 0;
        m_ifs.read(reinterpret_cast<char*>(&dim), sizeof(uint32_t));
        if (m_ifs.eof()) {
            m_vec.clear();
            return m_vec;
        }
        ABORT_IF_NE(m_vec.size(), dim);
        m_ifs.read(reinterpret_cast<char*>(m_vec.data()), sizeof(float) * dim);
        return m_vec;
    }

  private:
    std::ifstream m_ifs;
    std::vector<float> m_vec;
};

class bvecs_iterator {
  public:
    bvecs_iterator(const std::string& path, size_t dimension) : m_ifs(make_ifstream(path)), m_vec(dimension) {}

    const std::vector<float>& next() {
        uint32_t dim = 0;
        m_ifs.read(reinterpret_cast<char*>(&dim), sizeof(uint32_t));
        if (m_ifs.eof()) {
            m_vec.clear();
            return m_vec;
        }
        ABORT_IF_NE(m_vec.size(), dim);
        for (uint32_t i = 0; i < dim; ++i) {
            uint8_t v;
            m_ifs.read(reinterpret_cast<char*>(&v), sizeof(uint8_t));
            m_vec[i] = static_cast<float>(v);
        }
        return m_vec;
    }

  private:
    std::ifstream m_ifs;
    std::vector<float> m_vec;
};

class libsvm_iterator {
  public:
    libsvm_iterator(const std::string& path, size_t dimension) : m_ifs(make_ifstream(path)), m_vec(dimension) {}

    const std::vector<float>& next() {
        if (!std::getline(m_ifs, m_row)) {
            m_vec.clear();
            return m_vec;
        }

        std::fill(m_vec.begin(), m_vec.end(), 0.0);

        std::istringstream iss(m_row);
        iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');  // ignore the label

        for (std::string taken; iss >> taken;) {
            const auto pos = taken.find(':');
            ABORT_IF_EQ(pos, std::string::npos);

            const uint32_t id = static_cast<uint32_t>(std::stoul(taken.substr(0, pos))) - 1;
            const float weight = std::stof(taken.substr(pos + 1));

            ABORT_IF_LE(m_vec.size(), id);
            m_vec[id] = weight;
        }

        return m_vec;
    }

  private:
    std::string m_row;
    std::ifstream m_ifs;
    std::vector<float> m_vec;
};

class fasttext_iterator {
  public:
    fasttext_iterator(const std::string& path, size_t dimension) : m_ifs(make_ifstream(path)), m_vec(dimension) {}

    const std::vector<float>& next() {
        if (!std::getline(m_ifs, m_row)) {
            m_vec.clear();
            return m_vec;
        }

        std::fill(m_vec.begin(), m_vec.end(), 0.0);

        std::istringstream iss(m_row);
        for (std::string taken; iss >> taken;) {
            const auto pos = taken.find(':');
            ABORT_IF_EQ(pos, std::string::npos);

            const uint32_t id = static_cast<uint32_t>(std::stoul(taken.substr(0, pos)));
            const float weight = std::stof(taken.substr(pos + 1));

            ABORT_IF_LE(m_vec.size(), id);
            m_vec[id] = weight;
        }

        return m_vec;
    }

  private:
    std::string m_row;
    std::ifstream m_ifs;
    std::vector<float> m_vec;
};

std::vector<float> gen_random_matrix(size_t dimension, size_t repetitions, size_t seed) {
    std::vector<float> random_matrix(dimension * repetitions);
    std::default_random_engine engine(seed);
    std::normal_distribution<float> dist(0.0, 1.0);
    for (size_t i = 0; i < random_matrix.size(); ++i) {
        random_matrix[i] = dist(engine);
    }
    return random_matrix;
}

std::string get_ext(std::string path) {
    size_t pos = path.find_last_of(".") + 1;
    return path.substr(pos, path.size() - pos);
}

bool simhash(const float* data_vec, const float* random_vec, size_t dimension) {
    return std::inner_product(data_vec, data_vec + dimension, random_vec, 0.0) >= 0;
}

template <typename Iterator>
int main_template(const cmd_line_parser::parser& p) {
    auto input_path = p.get<std::string>("input_path");
    auto output_path = p.get<std::string>("output_path");
    auto dimension = p.get<size_t>("dimension");
    auto seed = p.get<size_t>("seed", 114514);

    Iterator data_itr(input_path, dimension);
    auto ofs = make_ofstream(tfm::format("%s.simhash%d", output_path, REPETITIONS));

    const std::vector<float> random_matrix = gen_random_matrix(dimension, REPETITIONS, seed);

    size_t num_codes = 0;
    for (;; ++num_codes) {
        if (num_codes % 100000 == 0) {
            tfm::printfln("%d vecs processed...", num_codes);
        }

        const std::vector<float>& data_vec = data_itr.next();
        if (data_vec.empty()) {
            break;
        }
        ABORT_IF_NE(data_vec.size(), dimension);

        uint64_t code = 0;
        for (size_t i = 0; i < REPETITIONS; ++i) {
            const float* random_vec = random_matrix.data() + (i * dimension);
            if (simhash(data_vec.data(), random_vec, dimension)) {
                code |= (1ULL << i);
            }
        }
        ofs.write(reinterpret_cast<const char*>(&code), sizeof(uint64_t));
    }

    tfm::printfln("Completed! The total number of vectors: %d", num_codes);
    return 0;
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "input file path of vectors");
    p.add("output_path", "output file path of hashed codes");
    p.add("format", "format (fvecs|bvecs|libsvm|fasttext)");
    p.add("dimension", "dimension of data vectors");
    p.add("seed", "seed", "-s", false);
    ABORT_IF(!p.parse());

    const auto format = p.get<std::string>("format");

    if (format == "fvecs") {
        return main_template<fvecs_iterator>(p);
    } else if (format == "bvecs") {
        return main_template<bvecs_iterator>(p);
    } else if (format == "libsvm") {
        return main_template<libsvm_iterator>(p);
    } else if (format == "fasttext") {
        return main_template<fasttext_iterator>(p);
    }

    tfm::errorfln("invalid format");
    return 1;
}