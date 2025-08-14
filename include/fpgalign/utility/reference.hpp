// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <fpgalign/config.hpp>

namespace utility
{

void store(std::vector<std::vector<uint8_t>> const & reference, config const & config, size_t const id);

void load(std::vector<std::vector<uint8_t>> & reference, config const & config, size_t const id);

} // namespace utility
