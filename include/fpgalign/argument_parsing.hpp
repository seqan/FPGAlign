// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint> // for uint8_t
#include <string>  // for string
#include <vector>  // for vector

#include <fpgalign/config.hpp> // for config

enum class subcommand : uint8_t
{
    build,
    search
};

struct parse_result
{
    subcommand subcmd;
    config cfg;
};

parse_result parse_arguments(std::vector<std::string> command_line);
