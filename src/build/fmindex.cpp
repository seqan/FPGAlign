// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/io/sequence_file/input.hpp>

#include <fmindex-collection/fmindex/BiFMIndex.h>

#include <fpgalign/build/build.hpp>

namespace build
{

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void fmindex(config const & config)
{
    std::vector<std::vector<uint8_t>> const reference = [&config]()
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        auto const bin_pathss = parse_input(config);

        std::vector<std::vector<uint8_t>> result{};

        for (auto && bin_paths : bin_pathss)
        {
            for (auto && bin_path : bin_paths)
            {
                sequence_file_t fin{bin_path};
                for (auto && record : fin)
                {
                    result.push_back({});
                    std::ranges::copy(record.sequence()
                                          | std::views::transform(
                                              [](auto const & in)
                                              {
                                                  return seqan3::to_rank(in);
                                              }),
                                      std::back_inserter(result.back()));
                }
            }
        }

        return result;
    }();

    fmc::BiFMIndex<4> index{reference, /*samplingRate*/ 16, config.threads};

    {
        std::ofstream os{config.output_path.string() + ".fmindex", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }
}

} // namespace build
