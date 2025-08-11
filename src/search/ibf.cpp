// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <random>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <hibf/config.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

#include <fpgalign/search/search.hpp>
#include <threshold/threshold.hpp>

namespace search
{

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

static inline constexpr uint64_t adjust_seed(uint8_t const kmer_size) noexcept
{
    return 0x8F3F73B5CF1C9ADEULL >> (64u - 2u * kmer_size);
}

threshold::threshold get_thresholder(config const & config)
{
    size_t const first_sequence_size = [&]()
    {
        seqan3::sequence_file_input<dna4_traits> fin{config.query_path};
        auto & record = *fin.begin();
        return record.sequence().size();
    }();

    return {threshold::threshold_parameters{.window_size = config.window_size,
                                            .shape = seqan3::ungapped{config.kmer_size},
                                            .query_length = first_sequence_size,
                                            .errors = config.errors}};
}

using seqfile_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>>;
using record_t = typename seqfile_t::record_type;

std::vector<hit> ibf(config const & config, size_t & todo_bin_count)
{
    seqan::hibf::interleaved_bloom_filter ibf{};

    {
        std::ifstream os{config.input_path.string() + ".ibf", std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        iarchive(ibf);
    }
    todo_bin_count = ibf.bin_count();

    std::vector<record_t> records = [&]()
    {
        std::vector<record_t> result{};
        seqfile_t fin{config.query_path};
        std::ranges::move(fin, std::back_inserter(result));
        // Very fast, improves parallel processing when chunks of the query belong to the same bin.
        std::ranges::shuffle(result, std::mt19937_64{0u});
        return result;
    }();

    std::vector<hit> hits(records.size());

#pragma omp parallel num_threads(config.threads)
    {
        auto agent = ibf.membership_agent();
        threshold::threshold const thresholder = get_thresholder(config);
        auto minimiser_view = seqan3::views::minimiser_hash(seqan3::ungapped{config.kmer_size},
                                                            seqan3::window_size{config.window_size},
                                                            seqan3::seed{adjust_seed(config.kmer_size)});
        std::vector<uint64_t> hashes;

#pragma omp for
        for (size_t i = 0; i < records.size(); ++i)
        {
            auto & [id, seq] = records[i];
            auto view = seq | minimiser_view | std::views::common;
            hashes.clear();
            hashes.assign(view.begin(), view.end());

            auto & result = agent.membership_for(hashes, thresholder.get(hashes.size()));
            hits[i].id = std::move(id);
            std::ranges::copy(std::move(seq)
                                  | std::views::transform(
                                      [](auto const & in)
                                      {
                                          return seqan3::to_rank(in);
                                      }),
                              std::back_inserter(hits[i].seq));
            std::ranges::copy(result, std::back_inserter(hits[i].bins));
        }
    }

    for (auto & hit : hits)
    {
        std::cout << hit.id << ": ";
        for (auto chr : hit.seq)
            std::cout << static_cast<uint16_t>(chr);
        std::cout << "\n    ";
        for (auto bin : hit.bins)
            std::cout << bin << ',';
        std::cout << '\n';
    }

    return hits;
}

} // namespace search
