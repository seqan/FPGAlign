// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <random>

#include <seqan3/io/sequence_file/input.hpp>

#include <hibf/config.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

#include <fpgalign/contrib/minimiser_hash.hpp>
#include <fpgalign/search/search.hpp>
#include <threshold/threshold.hpp>

namespace search
{

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

threshold::threshold get_thresholder(config const & config, meta const & meta)
{
    size_t const first_sequence_size = [&]()
    {
        seqan3::sequence_file_input<dna4_traits> fin{config.query_path};
        auto & record = *fin.begin();
        return record.sequence().size();
    }();

    return {threshold::threshold_parameters{.window_size = meta.window_size,
                                            .shape = seqan3::ungapped{meta.kmer_size},
                                            .query_length = first_sequence_size,
                                            .errors = config.errors}};
}

using seqfile_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>>;
using record_t = typename seqfile_t::record_type;

std::vector<hit> ibf(config const & config, meta & meta)
{
    seqan::hibf::interleaved_bloom_filter ibf{};

    {
        std::ifstream os{config.input_path.string() + ".ibf", std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        iarchive(ibf);
    }
    assert(ibf.bin_count() == meta.number_of_bins);

    meta.queries = [&]()
    {
        std::vector<record_t> result{};
        seqfile_t fin{config.query_path};
        std::ranges::move(fin, std::back_inserter(result));
        // Very fast, improves parallel processing when chunks of the query belong to the same bin.
        std::ranges::shuffle(result, std::mt19937_64{0u});
        return result;
    }();

    std::vector<hit> hits(meta.queries.size());

#pragma omp parallel num_threads(config.threads)
    {
        auto agent = ibf.membership_agent();
        threshold::threshold const thresholder = get_thresholder(config, meta);
        auto minimiser_view =
            contrib::views::minimiser_hash({.kmer_size = meta.kmer_size, .window_size = meta.window_size});

        std::vector<uint64_t> hashes;

#pragma omp for
        for (size_t i = 0; i < meta.queries.size(); ++i)
        {
            auto & [id, seq] = meta.queries[i];
            auto view = seq | minimiser_view | std::views::common;
            hashes.clear();
            hashes.assign(view.begin(), view.end());

            auto & result = agent.membership_for(hashes, thresholder.get(hashes.size()));
            std::ranges::copy(result, std::back_inserter(hits[i].bins));
        }
    }

    return hits;
}

} // namespace search
