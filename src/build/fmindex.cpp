// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/format.h>

#include <seqan3/io/sequence_file/input.hpp>

#include <hibf/contrib/std/enumerate_view.hpp>

#include <fmindex-collection/fmindex/BiFMIndex.h>

#include <fpgalign/build/build.hpp>

namespace build
{

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

std::vector<std::vector<uint8_t>> read_reference(std::vector<std::string> const & bin_paths)
{
    std::vector<std::vector<uint8_t>> reference{};

    for (auto const & bin_path : bin_paths)
    {
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{bin_path};
        for (auto && record : fin)
        {
            reference.push_back({});
            std::ranges::copy(record.sequence()
                                  | std::views::transform(
                                      [](auto const & in)
                                      {
                                          return seqan3::to_rank(in) + 1u;
                                      }),
                              std::back_inserter(reference.back()));
        }
    }

    return reference;
}

void fmindex(config const & config)
{
    for (auto && [id, bin_paths] : seqan::stl::views::enumerate(parse_input(config)))
    {
        auto reference = read_reference(bin_paths);
        fmc::BiFMIndex<5> index{reference, /*samplingRate*/ 16, config.threads};

        {
            std::ofstream os{fmt::format("{}.{}.fmindex", config.output_path.c_str(), id), std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(index);
        }
    }
}

} // namespace build
