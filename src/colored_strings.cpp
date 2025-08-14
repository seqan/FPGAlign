// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstring> // for memcmp
#include <string>  // for string, basic_string

#include <fmt/color.h> // for color, fg, format

#include <sharg/detail/terminal.hpp> // for stderr_is_terminal

#include <fpgalign/colored_strings.hpp> // for colored_strings

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
