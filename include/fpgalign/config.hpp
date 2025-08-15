// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>    // for uint8_t, uint16_t, uint32_t
#include <filesystem> // for path

struct config
{
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};

    uint8_t hash_count{2u};
    double fpr{0.05};

    std::filesystem::path input_path{};
    std::filesystem::path output_path{};
    std::filesystem::path query_path{};
    uint8_t errors{0u};
    uint16_t threads{1u};
    size_t queue_capacity{1u};
};
