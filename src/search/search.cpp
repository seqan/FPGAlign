// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sam_file/output.hpp>

#include <fpgalign/search/search.hpp>

namespace search
{

auto format_as(wip_alignment wip)
{
    return fmt::format("(bin = {}, seqNo = {}, pos = {}, seq = {:d}, ref = {:d})",
                       wip.bin,
                       wip.sequence_number,
                       wip.position,
                       fmt::join(wip.seq, ""),
                       fmt::join(wip.ref, ""));
}

void do_alignment(config const & config, std::vector<wip_alignment> const & wips)
{
    seqan3::sam_file_output sam_out{config.output_path,
                                    seqan3::fields<seqan3::field::seq,
                                                   seqan3::field::id,
                                                   seqan3::field::ref_id,
                                                   seqan3::field::ref_offset,
                                                   seqan3::field::cigar,
                                                   //    seqan3::field::qual,
                                                   seqan3::field::mapq>{}};

    seqan3::configuration const align_config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
        | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};

    for (auto & wip : wips)
    {
        size_t const start = wip.position - static_cast<size_t>(wip.position != 0u);
        size_t const length = wip.seq.size();
        auto it = std::ranges::next(wip.ref.begin(), start, wip.ref.end());
        auto end = std::ranges::next(it, length + 1u, wip.ref.end());
        std::span ref_text{it, end};

        for (auto && alignment : seqan3::align_pairwise(std::tie(ref_text, wip.seq), align_config))
        {
            auto cigar = seqan3::cigar_from_alignment(alignment.alignment());
            size_t ref_offset = alignment.sequence1_begin_position() + 2 + start;
            size_t map_qual = 60u + alignment.score();

            sam_out.emplace_back(std::views::transform(wip.seq,
                                                       [](auto const in)
                                                       {
                                                           return seqan3::dna4{}.assign_rank(
                                                               static_cast<uint8_t>(in - 1u));
                                                       }),
                                 wip.id,
                                 fmt::format("{}", wip.sequence_number), // todo ref storage
                                 ref_offset,
                                 cigar,
                                 //  record.base_qualities(),
                                 map_qual);
        }
    }
}

void search(config const & config)
{
    size_t todo_bin_count{};
    std::vector<hit> hits = ibf(config, todo_bin_count);
    auto const res = fmindex(config, std::move(hits), todo_bin_count);
    do_alignment(config, res);
    // for (auto const & elem : res)
    // {
    //     fmt::print("{}\n", elem);
    // }
}

} // namespace search
