// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sam_file/output.hpp>

#include <fpgalign/search/search.hpp>

namespace search
{

using sam_out_t = seqan3::sam_file_output<seqan3::fields<seqan3::field::seq,
                                                         seqan3::field::id,
                                                         seqan3::field::ref_id,
                                                         seqan3::field::ref_offset,
                                                         seqan3::field::cigar,
                                                         //    seqan3::field::qual,
                                                         seqan3::field::mapq>>;

void task(meta & meta, std::span<wip_alignment> wip, sam_out_t & sam_out)
{
    static seqan3::configuration const align_config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
        | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};

    for (auto & [bin, sequence_number, position, idx] : wip)
    {
        auto & seq = meta.queries[idx].sequence();
        auto seq_view = std::views::transform(seq,
                                              [](seqan3::dna4 const in) -> uint8_t
                                              {
                                                  return in.to_rank() + 1u;
                                              });
        auto & seq_id = meta.queries[idx].id();
        auto & ref = meta.references[bin][sequence_number];
        auto & ref_id = meta.ref_ids[bin][sequence_number];

        size_t const start = position - static_cast<size_t>(position != 0u);
        size_t const length = seq.size();
        auto it = std::ranges::next(ref.begin(), start, ref.end());
        auto end = std::ranges::next(it, length + 1u, ref.end());
        std::span ref_text{it, end};

        for (auto && alignment : seqan3::align_pairwise(std::tie(ref_text, seq_view), align_config))
        {
            auto cigar = seqan3::cigar_from_alignment(alignment.alignment());
            size_t ref_offset = alignment.sequence1_begin_position() + 2 + start;
            size_t map_qual = 60u + alignment.score();

            sam_out.emplace_back(seq,
                                 seq_id,
                                 ref_id,
                                 ref_offset,
                                 cigar,
                                 //  record.base_qualities(),
                                 map_qual);
        }
    }
}

void do_alignment(config const & config, meta & meta, scq::slotted_cart_queue<wip_alignment> & alignment_queue)
{
    sam_out_t sam_out{config.output_path};

    while (true)
    {
        scq::cart_future<wip_alignment> cart = alignment_queue.dequeue();
        if (!cart.valid())
            return;

        task(meta, cart.get().second, sam_out);
    }
}

} // namespace search
