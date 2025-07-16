// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/io/sequence_file/input.hpp>

#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/search/search.h>

#include <fpgalign/search/search.hpp>

namespace search
{

void fmindex(config const & config, std::vector<hit> hits)
{
    fmc::BiFMIndex<4> index{};

    {
        std::ifstream os{config.input_path.string() + ".fmindex", std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        iarchive(index);
    }

#pragma omp parallel for num_threads(config.threads)
    for (size_t i = 0; i < hits.size(); ++i)
    {
        auto & [id, seq, bins] = hits[i];
        auto callback = [&](auto cursor, size_t)
        {
            for (auto j : cursor)
            {
                auto [entry, offset] = index.locate(j);
                auto [seqId, pos] = entry;
#pragma omp critical
                {
                    std::cout << '[' << id << "] found hit in seqNo " << seqId << " Pos " << pos + offset << '\n';
                }
            }
        };
        fmc::search<true>(index, seq, config.errors, callback);
        (void)bins;
    }
}

} // namespace search
