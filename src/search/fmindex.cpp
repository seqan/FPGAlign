// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/search/search.h>

#include <fpgalign/contrib/slotted_cart_queue.hpp>
#include <fpgalign/search/search.hpp>
#include <fpgalign/utility/fmindex.hpp>

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
                        alignment_queue.enqueue(scq::slot_id{0u},
                                                alignment_info{.bin = slot.value,
                                                               .sequence_number = seqId,
                                                               .position = pos + offset,
                                                               .idx = idx});
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
