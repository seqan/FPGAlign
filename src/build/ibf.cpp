// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/io/sequence_file/input.hpp>

#include <hibf/config.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

#include <fpgalign/build/build.hpp>
#include <fpgalign/colored_strings.hpp>
#include <fpgalign/contrib/minimiser_hash.hpp>
#include <fpgalign/utility/ibf.hpp>

namespace build
{

void ibf(config const & config, meta & meta)
{
    meta.kmer_size = config.kmer_size;
    meta.window_size = config.window_size;

    auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        auto minimiser_view = contrib::views::minimiser_hash({.kmer_size = config.kmer_size, //
                                                              .window_size = config.window_size});

        for (auto && bin_path : meta.bin_paths[user_bin_id])
        {
            seqfile_t fin{bin_path};
            for (auto && record : fin)
            {
                if (size_t const record_size = record.sequence().size(); record_size < config.window_size)
                {
#pragma omp critical
                    {
                        std::cerr << colored_strings::cerr::warning << "File " << std::quoted(bin_path)
                                  << " contains a sequence of length " << record_size << " (ID=" << record.id()
                                  << "). This is shorter than the window size (" << config.window_size
                                  << ") and will result in no k-mers being generated for this sequence. A user bin "
                                     "without k-mers will result in an error.\n";
                    }
                }
                std::ranges::copy(record.sequence() | minimiser_view, it);
            }
        }
    };

    seqan::hibf::config ibf_config{.input_fn = get_user_bin_data,
                                   .number_of_user_bins = meta.bin_paths.size(),
                                   .number_of_hash_functions = config.hash_count,
                                   .maximum_fpr = config.fpr,
                                   .threads = config.threads};

    seqan::hibf::interleaved_bloom_filter ibf{ibf_config};

    utility::store(ibf, config);
}

} // namespace build
