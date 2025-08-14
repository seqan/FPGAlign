// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __move, __shuffle, move, shuffle
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t
#include <filesystem> // for path
#include <iterator>   // for back_insert_iterator, back_inserter, operator==
#include <random>     // for mt19937_64
#include <ranges>     // for common_view, operator|, __fn, common, views
#include <tuple>      // for get
#include <vector>     // for vector

#include <seqan3/contrib/std/detail/adaptor_base.hpp> // for operator|
#include <seqan3/io/detail/misc.hpp>                  // for set_format
#include <seqan3/io/record.hpp>                       // for fields
#include <seqan3/io/sequence_file/record.hpp>         // for sequence_record
#include <seqan3/search/kmer_index/shape.hpp>         // for shape, ungapped

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter

#include <fpgalign/config.hpp>                     // for config
#include <fpgalign/contrib/minimiser_hash.hpp>     // for minimiser_hash, operator==, minimiser_hash_fn
#include <fpgalign/contrib/slotted_cart_queue.hpp> // for slotted_cart_queue, assert, slot_id
#include <fpgalign/meta.hpp>                       // for meta, seqfile_t, record_t
#include <fpgalign/search/search.hpp>              // for ibf
#include <fpgalign/utility/ibf.hpp>                // for load
#include <threshold/threshold.hpp>                 // for threshold
#include <threshold/threshold_parameters.hpp>      // for threshold_parameters

namespace search
{

threshold::threshold get_thresholder(config const & config, meta const & meta)
{
    size_t const first_sequence_size = [&]()
    {
        seqfile_t fin{config.query_path};
        auto & record = *fin.begin();
        return record.sequence().size();
    }();

    return {threshold::threshold_parameters{.window_size = meta.window_size,
                                            .shape = seqan3::ungapped{meta.kmer_size},
                                            .query_length = first_sequence_size,
                                            .errors = config.errors}};
}

void ibf(config const & config, meta & meta, scq::slotted_cart_queue<size_t> & filter_queue)
{
    seqan::hibf::interleaved_bloom_filter ibf{};
    utility::load(ibf, config);

    assert(ibf.bin_count() == meta.number_of_bins);

    meta.queries = [&]()
    {
        std::vector<record_t> result{};
        seqfile_t fin{config.query_path};
        std::ranges::move(fin, std::back_inserter(result));
        // Very fast, improves parallel processing when chunks of the query belong to the same bin.
        std::ranges::shuffle(result, std::mt19937_64{0u});
        return result;
    }();

#pragma omp parallel num_threads(config.threads)
    {
        auto agent = ibf.membership_agent();
        threshold::threshold const thresholder = get_thresholder(config, meta);
        auto minimiser_view = contrib::views::minimiser_hash({.kmer_size = meta.kmer_size, //
                                                              .window_size = meta.window_size});

        std::vector<uint64_t> hashes;

#pragma omp for
        for (size_t i = 0; i < meta.queries.size(); ++i)
        {
            auto & [id, seq] = meta.queries[i];
            auto view = seq | minimiser_view | std::views::common;
            hashes.clear();
            hashes.assign(view.begin(), view.end());

            auto & result = agent.membership_for(hashes, thresholder.get(hashes.size()));
            for (size_t bin : result)
            {
                filter_queue.enqueue(scq::slot_id{bin}, i);
            }
        }
    }

    filter_queue.close();
}

} // namespace search
