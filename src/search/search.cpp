// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/format.h>

#include <fpgalign/search/search.hpp>

namespace search
{

void search(config const & config)
{
    meta meta{};
    {
        std::ifstream is{fmt::format("{}.meta", config.input_path.c_str()), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(meta);
    }

    meta.references.resize(meta.number_of_bins);
    for (size_t i = 0; i < meta.number_of_bins; ++i)
    {

        std::ifstream is{fmt::format("{}.{}.ref", config.input_path.c_str(), i), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(meta.references[i]);
    }

    std::vector<hit> hits = ibf(config, meta);
    auto const res = fmindex(config, meta, std::move(hits));
    do_alignment(config, meta, res);
}

} // namespace search
