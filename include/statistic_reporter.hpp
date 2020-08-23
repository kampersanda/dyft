#pragma once

#include <iomanip>
#include <map>
#include <sstream>
#include <variant>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant.hpp>
#include <boost/variant/get.hpp>

#include "tinyformat/tinyformat.h"

namespace dyft {

class statistic_reporter {
  public:
    using value_type = boost::variant<int,  // 0
                                      float,  // 1
                                      double,  // 2
                                      uint32_t,  // 3
                                      uint64_t,  // 4
                                      size_t,  // 5
                                      std::string,  // 6
                                      std::vector<int>,  // 7
                                      std::vector<float>,  // 8
                                      std::vector<double>,  // 9
                                      std::vector<uint32_t>,  // 10
                                      std::vector<uint64_t>,  // 11
                                      std::vector<size_t>  // 12
                                      >;
    using row_type = std::map<std::string, value_type>;
    using table_type = std::vector<row_type>;

    template <int N, typename... Ts>
    using NthTypeOf = typename std::tuple_element<N, std::tuple<Ts...>>::type;

    template <int N, typename... Ts>
    static auto& access(boost::variant<Ts...>& v) {
        using target = NthTypeOf<N, Ts...>;
        return boost::get<target>(v);
    }

    template <int N, typename... Ts>
    static auto& access(const boost::variant<Ts...>& v) {
        using target = NthTypeOf<N, Ts...>;
        return boost::get<target>(v);
    }

  private:
    boost::posix_time::ptime m_date;
    std::map<std::string, value_type> m_tags;
    std::map<std::string, table_type> m_tables;

  public:
    statistic_reporter() : m_date(boost::posix_time::second_clock::local_time()) {}

    static statistic_reporter& get_instance() {
        static statistic_reporter instance;
        return instance;
    }

    void tag(const std::string& key, const value_type& val) {
        m_tags[key] = val;
    }
    void append(const std::string& key, const row_type& row) {
        auto itr = m_tables.find(key);
        if (itr == m_tables.end()) {
            m_tables[key] = table_type{row};
        } else {
            itr->second.push_back(row);
        }
    }

    template <class T>
    static boost::property_tree::ptree vec_to_ptree(const std::vector<T>& vec) {
        boost::property_tree::ptree pt;
        for (auto v : vec) {
            boost::property_tree::ptree ct;
            ct.put("", v);
            pt.push_back(std::make_pair("", ct));
        }
        return pt;
    }

    template <size_t N = 0>
    static void update_ptree(boost::property_tree::ptree& pt, const std::string& key, const value_type& val) {
        if constexpr (N <= 6) {
            if (val.which() == N) {
                pt.put(key, access<N>(val));
            } else {
                update_ptree<N + 1>(pt, key, val);
            }
        } else if constexpr (N <= 12) {
            if (val.which() == N) {
                pt.add_child(key, vec_to_ptree(access<N>(val)));
            } else {
                update_ptree<N + 1>(pt, key, val);
            }
        }
    }

    boost::property_tree::ptree make_ptree() const {
        boost::property_tree::ptree root;
        root.put("date", boost::posix_time::to_iso_extended_string(m_date));
        for (const auto& v : m_tags) {
            update_ptree<0>(root, tfm::format("tags.%s", v.first), v.second);
        }
        for (const auto& t : m_tables) {
            boost::property_tree::ptree c;
            for (const auto& r : t.second) {
                boost::property_tree::ptree cc;
                for (const auto& v : r) {
                    update_ptree<0>(cc, tfm::format("%s", v.first), v.second);
                }
                c.push_back(std::make_pair("", cc));
            }
            root.add_child(t.first, c);
        }
        return root;
    }

    void save_json(std::string path) const {
        boost::property_tree::write_json(path, make_ptree());
    }
};

#define STATISTIC_TAG(key, val) statistic_reporter::get_instance().tag(std::string(key), val)
#define STATISTIC_APPEND(key, ...) \
    statistic_reporter::get_instance().append(std::string(key), statistic_reporter::row_type(__VA_ARGS__))
#define STATISTIC_SAVE(path) statistic_reporter::get_instance().save_json(path)

}  // namespace dyft
