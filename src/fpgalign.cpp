// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>

#include <fpgalign/argument_parsing.hpp>
#include <fpgalign/build/build.hpp>
#include <fpgalign/colored_strings.hpp>
#include <fpgalign/config.hpp>
#include <fpgalign/search/search.hpp>

int main(int argc, char ** argv)
{
    try
    {
        parse_result const result = parse_arguments({argv, argv + argc});

        if (result.subcmd == subcommand::build)
            build::build(result.cfg);
        if (result.subcmd == subcommand::search)
            search::search(result.cfg);
    }
    catch (std::exception const & ext)
    {
        std::cerr << colored_strings::cerr::error << ext.what() << '\n';
        std::exit(-1);
    }

    return 0;
}
