// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>

#include <fmt/format.h>

#include <cereal/archives/binary.hpp>

#include <fpgalign/utility/ibf.hpp>

namespace utility
{

void store(seqan::hibf::interleaved_bloom_filter const & ibf, config const & config)
{
    std::ofstream os{fmt::format("{}.ibf", config.input_path.c_str()), std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

void load(seqan::hibf::interleaved_bloom_filter & ibf, config const & config)
{
    std::ifstream is{fmt::format("{}.ibf", config.input_path.string()), std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}

} // namespace utility
