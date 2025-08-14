// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <fmindex-collection/fmindex/BiFMIndex.h>

#include <fpgalign/config.hpp>

namespace utility
{

void store(fmc::BiFMIndex<5> const & index, config const & config, size_t const id);

void load(fmc::BiFMIndex<5> & index, config const & config, size_t const id);

} // namespace utility
