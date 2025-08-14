// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert> // for assert
#include <fstream> // for char_traits, basic_istream, basic_ifstream, getline, operator>>, ifstream
#include <sstream> // for basic_istringstream
#include <string>  // for basic_string, string
#include <utility> // for move
#include <vector>  // for vector

#include <fpgalign/build/build.hpp>  // for fmindex, ibf, build, parse_input
#include <fpgalign/config.hpp>       // for config
#include <fpgalign/meta.hpp>         // for meta
#include <fpgalign/utility/meta.hpp> // for store

namespace build
{

std::vector<std::vector<std::string>> parse_input(config const & config)
{
    std::ifstream input_file{config.input_path};

    std::vector<std::vector<std::string>> result;
    std::string line;
    std::string token;

    while (std::getline(input_file, line))
    {
        if (line.empty())
            continue;

        std::vector<std::string> tokens;
        std::istringstream line_stream{line};

        while (line_stream >> token)
            tokens.push_back(token);

        result.push_back(std::move(tokens));
    }

    return result;
}

void build(config const & config)
{
    meta meta{};
    meta.bin_paths = parse_input(config);
    meta.number_of_bins = meta.bin_paths.size();
    build::ibf(config, meta);
    assert(meta.kmer_size == config.kmer_size);
    assert(meta.window_size == config.window_size);
    build::fmindex(config, meta);

    utility::store(meta, config);
}

} // namespace build
