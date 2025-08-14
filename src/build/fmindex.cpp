// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __copy, copy
#include <cstddef>    // for size_t
#include <cstdint>    // for uint8_t
#include <filesystem> // for path
#include <iterator>   // for back_insert_iterator, back_inserter
#include <ranges>     // for transform_view, views, operator|, __pipeable, __fn, operat...
#include <string>     // for basic_string
#include <tuple>      // for get
#include <vector>     // for vector

#include <seqan3/alphabet/concept.hpp> // for to_rank
#include <seqan3/io/detail/misc.hpp>   // for set_format
#include <seqan3/io/record.hpp>        // for fields, field

#include <fmindex-collection/fmindex/BiFMIndex.h> // for BiFMIndex

#include <fpgalign/build/build.hpp>       // for fmindex
#include <fpgalign/config.hpp>            // for config
#include <fpgalign/meta.hpp>              // for meta, seqfile_t
#include <fpgalign/utility/fmindex.hpp>   // for store
#include <fpgalign/utility/reference.hpp> // for store
#include <std/detail/adaptor_base.hpp>    // for operator|

namespace build
{

void read_reference_into(std::vector<std::vector<uint8_t>> & reference, meta & meta, size_t const i)
{
    reference.clear();

    for (auto const & bin_path : meta.bin_paths[i])
    {
        seqfile_t fin{bin_path};

        for (auto && record : fin)
        {
            meta.ref_ids[i].push_back(record.id());
            reference.push_back({});
            std::ranges::copy(record.sequence()
                                  | std::views::transform(
                                      [](auto const & in)
                                      {
                                          return seqan3::to_rank(in) + 1u;
                                      }),
                              std::back_inserter(reference.back()));
        }
    }
}

void fmindex(config const & config, meta & meta)
{
    meta.ref_ids.resize(meta.number_of_bins);

#pragma omp parallel num_threads(config.threads)
    {
        std::vector<std::vector<uint8_t>> reference;

#pragma omp for
        for (size_t i = 0; i < meta.number_of_bins; ++i)
        {
            read_reference_into(reference, meta, i);

            fmc::BiFMIndex<5> index{reference, /*samplingRate*/ 16, /*threads*/ 1u};

            utility::store(index, config, i);
            utility::store(reference, config, i);
        }
    }
}

} // namespace build
