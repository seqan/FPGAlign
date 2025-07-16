// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <hibf/config.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

#include <fpgalign/build/build.hpp>
#include <fpgalign/colored_strings.hpp>

namespace build
{

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

static inline constexpr uint64_t adjust_seed(uint8_t const kmer_size) noexcept
{
    return 0x8F3F73B5CF1C9ADEULL >> (64u - 2u * kmer_size);
}

void ibf(config const & config)
{
    auto const bin_paths = parse_input(config);

    auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        auto minimiser_view = seqan3::views::minimiser_hash(seqan3::ungapped{config.kmer_size},
                                                            seqan3::window_size{config.window_size},
                                                            seqan3::seed{adjust_seed(config.kmer_size)});

        for (auto && bin_path : bin_paths[user_bin_id])
        {
            sequence_file_t fin{bin_path};
            for (auto && record : fin)
            {
                if (size_t const record_size = record.sequence().size(); record_size < config.window_size)
                {
#pragma omp critical
                    {
                        std::cerr << colored_strings::cerr::warning << "File " << std::quoted(bin_path)
                                  << " contains a sequence of length " << record_size
                                  << ". This is shorter than the window size (" << config.window_size
                                  << ") and will result in no k-mers being generated for this sequence. A user bin "
                                     "without k-mers will result in an error.\n";
                    }
                }
                std::ranges::copy(record.sequence() | minimiser_view, it);
            }
        }
    };

    seqan::hibf::config ibf_config{.input_fn = get_user_bin_data,
                                   .number_of_user_bins = bin_paths.size(),
                                   .number_of_hash_functions = config.hash_count,
                                   .maximum_fpr = config.fpr,
                                   .threads = config.threads};

    seqan::hibf::interleaved_bloom_filter ibf{ibf_config};

    {
        std::ofstream os{config.output_path.string() + ".ibf", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(ibf);
    }
}

} // namespace build
