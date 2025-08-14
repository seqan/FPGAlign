// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>

#include <fmt/format.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <fpgalign/utility/reference.hpp>

namespace utility
{

void store(std::vector<std::vector<uint8_t>> const & reference, config const & config, size_t const id)
{
    std::ofstream os{fmt::format("{}.{}.ref", config.output_path.c_str(), id), std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(reference);
}

void load(std::vector<std::vector<uint8_t>> & reference, config const & config, size_t const id)
{
    std::ifstream is{fmt::format("{}.{}.ref", config.input_path.c_str(), id), std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(reference);
}

} // namespace utility
