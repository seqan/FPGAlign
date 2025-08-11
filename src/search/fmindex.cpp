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

fmc::BiFMIndex<4> load_index(config const & config, size_t const id)
{
    fmc::BiFMIndex<4> index{};

    {
        std::ifstream os{fmt::format("{}.{}.fmindex", config.input_path.c_str(), id), std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        iarchive(index);
    }

    return index;
}

void fmindex(config const & config, std::vector<hit> hits, size_t const todo_bin_count)
{
    // todo bin count
    // todo capacity
    // each slot = 1 bin
    // a cart is full if it has 5 elements (hits)
    scq::slotted_cart_queue<size_t> queue{{.slots = todo_bin_count, .carts = todo_bin_count, .capacity = 5}};
    size_t thread_id{};

    auto get_thread = [&]()
    {
        return std::jthread(
            [&, thread_id = thread_id++]()
            {
                while (true)
                {
                    scq::cart_future<size_t> cart = queue.dequeue();
                    if (!cart.valid())
                        return;
                    auto [slot, span] = cart.get();
                    auto index = load_index(config, slot.value);
                    for (auto idx : span)
                    {
                        auto & [id, seq, bins] = hits[idx];

                        auto callback = [&](auto cursor, size_t)
                        {
                            for (auto j : cursor)
                            {
                                auto [entry, offset] = index.locate(j);
                                auto [seqId, pos] = entry;
                                {
                                    fmt::print("[{}][{}] found hit in bin {} in seqNo {} at Pos {}\n",
                                               thread_id,
                                               id,
                                               slot.value,
                                               seqId,
                                               pos + offset);
                                }
                            }
                        };

                        fmc::search<true>(index, seq, config.errors, callback);
                    }
                }
            });
    };

    std::vector<std::jthread> worker(config.threads);
    std::ranges::generate(worker, get_thread);

    for (auto && [idx, hit] : seqan::stl::views::enumerate(hits))
        for (auto bin : hit.bins)
            queue.enqueue(scq::slot_id{bin}, idx);

    queue.close();
}

} // namespace search
