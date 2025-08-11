// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <sharg/parser.hpp>

#include <fpgalign/argument_parsing.hpp>

namespace build
{

config parse_arguments(sharg::parser & parser)
{
    config config{};

    parser.add_subsection("General options");
    parser.add_option(config.input_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "A file containing file names",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(config.output_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Output path",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}}); // .ibf and .fmindex
    parser.add_option(config.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use."}); // positive_integer_validator

    parser.add_subsection("k-mer options");
    parser.add_option(config.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(config.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The window size.",
                                    .default_message = "k-mer size"});

    parser.add_subsection("IBF options");
    parser.add_option(config.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The false positive rate.",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
    parser.add_option(config.hash_count,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "The number of hash functions to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 5}});

    parser.parse();

    if (parser.is_option_set("kmer") && !parser.is_option_set("window"))
        config.window_size = config.kmer_size;
    else if (config.window_size < config.kmer_size)
        throw sharg::validation_error{"k-mer size cannot be smaller than window size!"};

    return config;
}

} // namespace build

namespace search
{

config parse_arguments(sharg::parser & parser)
{
    config config{};

    parser.add_subsection("General options");
    parser.add_option(config.input_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "Prefix",
                                    .required = true}); // .ibf and .fmindex
    parser.add_option(config.query_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query",
                                    .description = "Query path",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(config.output_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Output path",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});
    parser.add_option(config.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use."}); // positive_integer_validator

    parser.add_subsection("k-mer options");
    parser.add_option(config.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(config.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The window size.",
                                    .default_message = "k-mer size"});
    parser.add_option(config.errors,
                      sharg::config{.short_id = '\0',
                                    .long_id = "errors",
                                    .description = "errors.",
                                    .validator = sharg::arithmetic_range_validator{0, 5}});

    parser.parse();

    if (parser.is_option_set("kmer") && !parser.is_option_set("window"))
        config.window_size = config.kmer_size;
    else if (config.window_size < config.kmer_size)
        throw sharg::validation_error{"k-mer size cannot be smaller than window size!"};

    return config;
}

} // namespace search

parse_result parse_arguments(std::vector<std::string> command_line)
{
    sharg::parser parser{"FPGAlign", std::move(command_line)};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";
    parser.add_subcommands({"build", "search"});

    parser.parse();

    sharg::parser & sub_parser = parser.get_sub_parser();
    parse_result result{};

    if (sub_parser.info.app_name == std::string_view{"FPGAlign-build"})
    {
        result.subcmd = subcommand::build;
        result.cfg = build::parse_arguments(sub_parser);
    }
    if (sub_parser.info.app_name == std::string_view{"FPGAlign-search"})
    {
        result.subcmd = subcommand::search;
        result.cfg = search::parse_arguments(sub_parser);
    }

    return result;
}
