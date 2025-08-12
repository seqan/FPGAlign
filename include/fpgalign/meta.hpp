// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <seqan3/io/sequence_file/input.hpp>

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using seqfile_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>>;
using record_t = typename seqfile_t::record_type;

struct meta
{
    uint8_t kmer_size{};
    uint32_t window_size{};
    size_t number_of_bins{};
    std::vector<std::vector<std::string>> bin_paths;
    std::vector<std::vector<std::string>> ref_ids;
    std::vector<std::vector<std::vector<uint8_t>>> references;
    std::vector<record_t> queries;

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(kmer_size);
        archive(window_size);
        archive(number_of_bins);
        archive(ref_ids);
    }
};
