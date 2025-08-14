// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <fpgalign/config.hpp>
#include <fpgalign/contrib/slotted_cart_queue.hpp>
#include <fpgalign/meta.hpp>

namespace search
{

struct alignment_info
{
    size_t bin;
    size_t sequence_number;
    size_t position;
    size_t idx;
};

void search(config const & config);
void ibf(config const & config, meta & meta, scq::slotted_cart_queue<size_t> & filter_queue);
void fmindex(config const & config,
             meta & meta,
             scq::slotted_cart_queue<size_t> & filter_queue,
             scq::slotted_cart_queue<alignment_info> & alignment_queue);
void do_alignment(config const & config, meta & meta, scq::slotted_cart_queue<alignment_info> & alignment_queue);

} // namespace search
