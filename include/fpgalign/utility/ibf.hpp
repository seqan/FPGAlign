// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <hibf/interleaved_bloom_filter.hpp>

#include <fpgalign/config.hpp>

namespace utility
{

void store(seqan::hibf::interleaved_bloom_filter const & ibf, config const & config);

void load(seqan::hibf::interleaved_bloom_filter & ibf, config const & config);

} // namespace utility
