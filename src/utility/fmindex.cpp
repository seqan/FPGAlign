// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>

#include <fmt/format.h>

#include <cereal/archives/binary.hpp>

#include <fpgalign/utility/fmindex.hpp>

namespace utility
{

void store(fmc::BiFMIndex<5> const & index, config const & config, size_t const id)
{
    std::ofstream os{fmt::format("{}.{}.fmindex", config.output_path.c_str(), id), std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

void load(fmc::BiFMIndex<5> & index, config const & config, size_t const id)
{
    std::ifstream is{fmt::format("{}.{}.fmindex", config.input_path.c_str(), id), std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(index);
}

} // namespace utility
