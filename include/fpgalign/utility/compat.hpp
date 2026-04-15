// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <fmindex-collection/search/SearchNg26.h>

namespace fmc::search_ng26
{

template <bool Edit = true, typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const & index,
            query_t && query,
            size_t maxErrors,
            delegate_t && delegate,
            size_t n = std::numeric_limits<size_t>::max())
{
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto
    {
        auto const & search_scheme = getCachedSearchScheme<Edit>(0, maxErrors, /*.shortLen=*/(length == 2));
        auto const & partition = getCachedPartition(search_scheme[0].pi.size(), length);
        return std::tie(search_scheme, partition);
    };
    size_t ct{};
    auto const & [search_scheme, partition] = selectSearchScheme(query.size());
    search_impl<Edit>(index,
                      query,
                      search_scheme,
                      partition,
                      [&](auto cur, size_t e)
                      {
                          if (cur.count() + ct > n)
                          {
                              cur.len = n - ct;
                          }
                          ct += cur.count();
                          delegate(cur, e);
                          return ct == n;
                      });
}

} // namespace fmc::search_ng26
