// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstddef> // for size_t
#include <string>  // for basic_string
#include <thread>  // for jthread
#include <vector>  // for vector

#include <fpgalign/config.hpp>                     // for config
#include <fpgalign/contrib/slotted_cart_queue.hpp> // for slotted_cart_queue
#include <fpgalign/meta.hpp>                       // for meta
#include <fpgalign/search/search.hpp>              // for alignment_info, do_alignment, fmindex, ibf, search
#include <fpgalign/utility/meta.hpp>               // for load
#include <fpgalign/utility/reference.hpp>          // for load

namespace search
{

void search(config const & config)
{
    meta meta{};
    utility::load(meta, config);

    meta.references.resize(meta.number_of_bins);
    for (size_t i = 0; i < meta.number_of_bins; ++i)
        utility::load(meta.references[i], config, i);

    // todo capacity
    // each slot = 1 bin
    // a cart is full if it has capacity many elements (hits)
    scq::slotted_cart_queue<size_t> filter_queue{{.slots = meta.number_of_bins, //
                                                  .carts = meta.number_of_bins,
                                                  .capacity = 1u}};
    scq::slotted_cart_queue<alignment_info> alignment_queue{{.slots = 1u, //
                                                             .carts = 1u,
                                                             .capacity = 1u}};

    std::jthread ibf_thread(
        [&]()
        {
            ibf(config, meta, filter_queue);
        });
    std::jthread fmindex_thread(
        [&]()
        {
            fmindex(config, meta, filter_queue, alignment_queue);
        });

    do_alignment(config, meta, alignment_queue);
}

} // namespace search
