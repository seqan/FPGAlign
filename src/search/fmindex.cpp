// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <seqan3/io/sequence_file/input.hpp>

#include <hibf/contrib/std/enumerate_view.hpp>

#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/search/search.h>

#include <fpgalign/contrib/slotted_cart_queue.hpp>
#include <fpgalign/search/search.hpp>

namespace search
{

// template <typename Index>
// auto myReconstruct(Index const & index, size_t seqNbr) -> std::vector<uint8_t>
// {
//     //TODO: möglicher weise ist das identisch zu einfach index.C[1]
//     auto totalNumberOfSeq = index.bwt.rank(index.size(), 0) + index.C[0];
//     for (size_t i{0}; i < totalNumberOfSeq; ++i)
//     {
//         auto idx = std::get<0>(std::get<0>(index.locate(i)));
//         if (idx == seqNbr)
//         {
//             return reconstructText(index, i);
//         }
//     }
//     throw std::runtime_error{"unknown sequence number"};
// }

fmc::BiFMIndex<5> load_index(config const & config, size_t const id)
{
    fmc::BiFMIndex<5> index{};

    {
        std::ifstream os{fmt::format("{}.{}.fmindex", config.input_path.c_str(), id), std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        iarchive(index);
    }

    // {
    //     fmt::println(" === Reconstruct all ===");
    //     auto text = fmc::reconstructText(index);
    //     for (size_t i = 0; i < text.size(); ++i)
    //         fmt::print("### {} ###\n{}\n", i, fmt::join(text[i], ""));
    // }
    // {
    //     fmt::println(" === Reconstruct single ===");
    //     auto text = fmc::reconstructText(index, 0);
    //     fmt::print("### 0 ###\n{}\n", fmt::join(text, ""));
    //     text = fmc::reconstructText(index, 1);
    //     fmt::print("### 1 ###\n{}\n", fmt::join(text, ""));
    // }
    // {
    //     fmt::println(" === Reconstruct alternative ===");
    //     auto text = myReconstruct(index, 0);
    //     fmt::print("### 0 ###\n{}\n", fmt::join(text, ""));
    //     text = myReconstruct(index, 1);
    //     fmt::print("### 1 ###\n{}\n", fmt::join(text, ""));
    // }

    return index;
}

std::vector<wip_alignment> fmindex(config const & config, meta & meta, std::vector<hit> hits)
{
    // todo bin count
    // todo capacity
    // each slot = 1 bin
    // a cart is full if it has 5 elements (hits)
    alignment_vector res;
    {
        scq::slotted_cart_queue<size_t> queue{
            {.slots = meta.number_of_bins, .carts = meta.number_of_bins, .capacity = 5}};
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
                            auto callback = [&](auto cursor, size_t)
                            {
                                for (auto j : cursor)
                                {
                                    auto [entry, offset] = index.locate(j);
                                    auto [seqId, pos] = entry;
                                    res.emplace_back(wip_alignment{.bin = slot.value,
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
                });
        };

        std::vector<std::jthread> worker(config.threads);
        std::ranges::generate(worker, get_thread);

        for (auto && [idx, hit] : seqan::stl::views::enumerate(hits))
            for (auto bin : hit.bins)
                queue.enqueue(scq::slot_id{bin}, idx);

        queue.close();
    } // Wait for threads to finish

    return res.get();
}

} // namespace search
