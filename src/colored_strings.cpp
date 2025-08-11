// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmt/color.h>

#include <sharg/detail/terminal.hpp>

#include <fpgalign/colored_strings.hpp>

bool const colored_strings::cerr_is_terminal = sharg::detail::stderr_is_terminal();

std::string const colored_strings::cerr::error = []() -> std::string
{
    if (cerr_is_terminal)
        return fmt::format(fg(fmt::color::red), "[Error] ");
    else
        return "[Error] ";
}();

std::string const colored_strings::cerr::warning = []() -> std::string
{
    if (cerr_is_terminal)
        return fmt::format(fg(fmt::color::yellow), "[Warning] ");
    else
        return "[Warning] ";
}();
