// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstring>    // for memcmp
#include <filesystem> // for path
#include <fstream>    // for basic_ifstream, basic_ofstream, basic_ios, ios, ifstream, ofstream

#include <fmt/format.h> // for format

#include <cereal/archives/binary.hpp> // for BinaryInputArchive, BinaryOutputArchive
#include <cereal/types/string.hpp>    // IWYU pragma: keep
#include <cereal/types/vector.hpp>    // IWYU pragma: keep

#include <fpgalign/config.hpp>       // for config
#include <fpgalign/meta.hpp>         // for meta
#include <fpgalign/utility/meta.hpp> // for load, store

namespace utility
{

void store(meta const & meta, config const & config)
{
    std::ofstream os{fmt::format("{}.meta", config.output_path.c_str()), std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(meta);
}

void load(meta & meta, config const & config)
{
    std::ifstream is{fmt::format("{}.meta", config.input_path.c_str()), std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(meta);
}

} // namespace utility
