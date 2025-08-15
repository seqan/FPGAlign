// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstddef> // for size_t
#include <cstdint> // for uint8_t
#include <ranges>  // for transform_view, __fn, transform, views
#include <tuple>   // for get
#include <utility> // for get
#include <vector>  // for vector

#include <seqan3/alphabet/nucleotide/dna4.hpp> // for dna4
#include <seqan3/io/sequence_file/record.hpp>  // for sequence_record

#include <fmindex-collection/fmindex/BiFMIndex.h> // for BiFMIndex
#include <fmindex-collection/search/search.h>     // for search

#include <fpgalign/config.hpp>                     // for config
#include <fpgalign/contrib/slotted_cart_queue.hpp> // for slotted_cart_queue, cart_future, slot_id, span
#include <fpgalign/meta.hpp>                       // for meta
#include <fpgalign/search/search.hpp>              // for alignment_info, fmindex
#include <fpgalign/utility/fmindex.hpp>            // for load

namespace search
{

void fmindex(config const & config,
             meta & meta,
             scq::slotted_cart_queue<size_t> & filter_queue,
             scq::slotted_cart_queue<alignment_info> & alignment_queue)
{
#pragma omp parallel num_threads(config.threads)
    {
        while (true)
        {
            scq::cart_future<size_t> cart = filter_queue.dequeue();
            if (!cart.valid())
                break;
            auto [slot, span] = cart.get();
            fmc::BiFMIndex<5> index{};
            utility::load(index, config, slot.value);
            for (auto idx : span)
            {
                auto callback = [&](auto cursor, size_t)
                {
                    for (auto j : cursor)
                    {
                        auto [entry, offset] = index.locate(j);
                        auto [seqId, pos] = entry;
                        alignment_queue.enqueue(slot,
                                                alignment_info{.query_idx = idx,
                                                               .reference_number = seqId,
                                                               .reference_position = pos + offset});
                    }
                };

                auto seq_view = std::views::transform(meta.queries[idx].sequence(),
                                                      [](seqan3::dna4 const in) -> uint8_t
                                                      {
                                                          return in.to_rank() + 1u;
                                                      });

                fmc::search<true>(index, seq_view, config.errors, callback);
            }
        }
    }

    alignment_queue.close();
}

} // namespace search
