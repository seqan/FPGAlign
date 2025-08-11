// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

struct colored_strings
{
    static bool const cerr_is_terminal;

    struct cerr
    {
        static std::string const error;
        static std::string const warning;
    };
};
