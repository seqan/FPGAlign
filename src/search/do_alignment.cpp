// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for find_if
#include <charconv>   // for to_chars
#include <cstddef>    // for size_t
#include <cstdint>    // for uint8_t
#include <filesystem> // for path
#include <iterator>   // for __next, next
#include <ranges>     // for transform_view, __fn, tra...
#include <string>     // for basic_string
#include <tuple>      // for tuple, tuple_cat, tie
#include <utility>    // for pair
#include <vector>     // for vector

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>              // for cigar_from_alignment
#include <seqan3/alignment/configuration/align_config_edit.hpp>                    // for edit_scheme
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>         // for gap_cost_affine
#include <seqan3/alignment/configuration/align_config_method.hpp>                  // for method_global, free_end_g...
#include <seqan3/alignment/configuration/align_config_output.hpp>                  // for output_alignment, output_...
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>          // for scoring_scheme
#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>     // for operator==
#include <seqan3/alignment/matrix/detail/trace_iterator_base.hpp>                  // for operator!=
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_base.hpp> // for matrix_major_order
#include <seqan3/alignment/pairwise/align_pairwise.hpp>                            // for align_pairwise
#include <seqan3/alignment/pairwise/alignment_result.hpp>                          // for alignment_result
#include <seqan3/alphabet/alphabet_base.hpp>                                       // for operator==, operator<
#include <seqan3/alphabet/nucleotide/dna4.hpp>                                     // for dna4
#include <seqan3/contrib/std/chunk_view.hpp>                                       // for operator==
#include <seqan3/contrib/std/detail/adaptor_base.hpp>                              // for operator|
#include <seqan3/contrib/std/pair.hpp>                                             // for get
#include <seqan3/contrib/std/tuple.hpp>                                            // for get
#include <seqan3/contrib/std/zip_view.hpp>                                         // for operator==
#include <seqan3/core/add_enum_bitwise_operators.hpp>                              // for operator|, operator&, ope...
#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>              // for algorithm_result_generato...
#include <seqan3/core/algorithm/detail/algorithm_executor_blocking.hpp>            // for algorithm_executor_blocking
#include <seqan3/core/configuration/configuration.hpp>                             // for configuration, operator|
#include <seqan3/core/range/detail/adaptor_base.hpp>                               // for operator|
#include <seqan3/io/detail/misc.hpp>                                               // for set_format
#include <seqan3/io/record.hpp>                                                    // for field, fields
#include <seqan3/io/sam_file/output.hpp>                                           // for sam_file_output
#include <seqan3/io/sequence_file/record.hpp>                                      // for sequence_record

#include <fpgalign/config.hpp>                     // for config
#include <fpgalign/contrib/slotted_cart_queue.hpp> // for span, cart_future, slotte...
#include <fpgalign/meta.hpp>                       // for meta
#include <fpgalign/search/search.hpp>              // for alignment_info, do_alignment

namespace search
{

using sam_out_t = seqan3::sam_file_output<seqan3::fields<seqan3::field::seq,
                                                         seqan3::field::id,
                                                         seqan3::field::ref_id,
                                                         seqan3::field::ref_offset,
                                                         seqan3::field::cigar,
                                                         //    seqan3::field::qual,
                                                         seqan3::field::mapq>>;

void task(meta & meta, size_t const bin, std::span<alignment_info> alignment_infos, sam_out_t & sam_out)
{
    static seqan3::configuration const align_config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
        | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};

    for (auto [query_idx, reference_number, reference_position] : alignment_infos)
    {
        auto & seq = meta.queries[query_idx].sequence();
        auto seq_view = std::views::transform(seq,
                                              [](seqan3::dna4 const in) -> uint8_t
                                              {
                                                  return in.to_rank() + 1u;
                                              });
        auto & seq_id = meta.queries[query_idx].id();
        auto & ref = meta.references[bin][reference_number];
        auto & ref_id = meta.ref_ids[bin][reference_number];

        size_t const start = reference_position - static_cast<size_t>(reference_position != 0u);
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

void do_alignment(config const & config, meta & meta, scq::slotted_cart_queue<alignment_info> & alignment_queue)
{
    sam_out_t sam_out{config.output_path};

    while (true)
    {
        scq::cart_future<alignment_info> cart = alignment_queue.dequeue();
        if (!cart.valid())
            return;
        auto [bin, alignment_infos] = cart.get();
        task(meta, bin.value, alignment_infos, sam_out);
    }
}

} // namespace search
