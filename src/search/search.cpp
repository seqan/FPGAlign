// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fpgalign/search/search.hpp>

namespace search
{

void search(config const & config)
{
    std::vector<hit> hits = ibf(config);
    fmindex(config, std::move(hits));
}

} // namespace search
