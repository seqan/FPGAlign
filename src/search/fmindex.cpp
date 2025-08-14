// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/format.h>

#include <seqan3/io/sequence_file/input.hpp>

#include <hibf/contrib/std/enumerate_view.hpp>

#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/search/search.h>

#include <fpgalign/contrib/slotted_cart_queue.hpp>
#include <fpgalign/search/search.hpp>

namespace search
{

fmc::BiFMIndex<5> load_index(config const & config, size_t const id)
{
    fmc::BiFMIndex<5> index{};

    {
        std::ifstream os{fmt::format("{}.{}.fmindex", config.input_path.c_str(), id), std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        iarchive(index);
    }

    return index;
}

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
            auto index = load_index(config, slot.value);
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
