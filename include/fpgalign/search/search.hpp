// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

#include <fpgalign/config.hpp>

namespace search
{

struct hit
{
    std::string id;
    std::vector<uint8_t> seq;
    std::vector<uint64_t> bins;
};

void search(config const & config);
std::vector<hit> ibf(config const & config, size_t & todo_bin_count);
void fmindex(config const & config, std::vector<hit> hits, size_t const todo_bin_count);

} // namespace search
