// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>

#include <fmt/format.h>

#include <cereal/archives/binary.hpp>

#include <fpgalign/utility/meta.hpp>

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
