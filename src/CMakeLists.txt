# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

set (FPGAlign_SOURCE_FILES
        argument_parsing.cpp
        build/build.cpp
        build/fmindex.cpp
        build/ibf.cpp
        colored_strings.cpp
        search/ibf.cpp
        search/fmindex.cpp
        search/search.cpp
)

# An object library (without main) to be used in multiple targets.
# You can add more external include paths of other projects that are needed for your project.
add_library (FPGAlign_lib STATIC ${FPGAlign_SOURCE_FILES})

target_include_directories (FPGAlign_lib PUBLIC "${FPGAlign_SOURCE_DIR}/include")
target_link_libraries (FPGAlign_lib PUBLIC seqan3::seqan3 sharg::sharg seqan::hibf fmt::fmt fmindex-collection::fmindex-collection seqan::threshold)

target_compile_options (FPGAlign_lib PUBLIC "-pedantic" "-Wall" "-Wextra")

option (FPGAlign_WITH_WERROR "Report compiler warnings as errors." ON)
if (FPGAlign_WITH_WERROR)
    target_compile_options (FPGAlign_lib PUBLIC "-Werror")
    message (STATUS "Building tests with -Werror.")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
        target_compile_options (FPGAlign_lib PUBLIC "-Wno-interference-size")
    endif ()

    # Warn about failed return value optimization.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
        target_compile_options (FPGAlign_lib PUBLIC "-Wnrvo")
    endif ()
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # std::jthread is experimental in LLVM < 20
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 20)
        target_compile_definitions (FPGAlign_lib PUBLIC "_LIBCPP_ENABLE_EXPERIMENTAL")
    endif ()
endif ()

# Add the application.
add_executable (FPGAlign fpgalign.cpp)
target_link_libraries (FPGAlign PRIVATE FPGAlign_lib)
