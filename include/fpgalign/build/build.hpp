// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

#include <fpgalign/config.hpp>
#include <fpgalign/meta.hpp>

namespace build
{

std::vector<std::vector<std::string>> parse_input(config const & config);

void build(config const & config);
void ibf(config const & config, meta & meta);
void fmindex(config const & config, meta & meta);

} // namespace build
