// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/format.h>

#include <seqan3/io/sequence_file/input.hpp>

#include <fmindex-collection/fmindex/BiFMIndex.h>

#include <fpgalign/build/build.hpp>

namespace build
{

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void read_reference_into(std::vector<std::vector<uint8_t>> & reference, meta & meta, size_t const i)
{
    reference.clear();

    for (auto const & bin_path : meta.bin_paths[i])
    {
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> fin{bin_path};

        for (auto && record : fin)
        {
            meta.ref_ids[i].push_back(record.id());
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
}

void fmindex(config const & config, meta & meta)
{
    meta.ref_ids.resize(meta.number_of_bins);

#pragma omp parallel num_threads(config.threads)
    {
        std::vector<std::vector<uint8_t>> reference;

#pragma omp for
        for (size_t i = 0; i < meta.number_of_bins; ++i)
        {
            read_reference_into(reference, meta, i);

            fmc::BiFMIndex<5> index{reference, /*samplingRate*/ 16, config.threads};

            {
                std::ofstream os{fmt::format("{}.{}.fmindex", config.output_path.c_str(), i), std::ios::binary};
                cereal::BinaryOutputArchive oarchive{os};
                oarchive(index);
            }

            {
                std::ofstream os{fmt::format("{}.{}.ref", config.output_path.c_str(), i), std::ios::binary};
                cereal::BinaryOutputArchive oarchive{os};
                oarchive(reference);
            }
        }
    }
}

} // namespace build
