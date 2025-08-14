// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <fpgalign/config.hpp>
#include <fpgalign/meta.hpp>

namespace utility
{

void store(meta const & meta, config const & config);

void load(meta & meta, config const & config);

} // namespace utility
