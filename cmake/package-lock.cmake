# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# CPM Package Lock (https://github.com/cpm-cmake/CPM.cmake)
# This file should be committed to version control

# cmake-format: off

# hibf
set (HIBF_VERSION e60abe791268be4e384740c52d29edc18606ea5f CACHE STRING "")
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${HIBF_VERSION} # main
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# sharg
set (SHARG_VERSION 47d4acc17a1e241ab2de9bb99d2d42a21f0d02ea CACHE STRING "")
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${SHARG_VERSION} # main
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "SHARG_NO_TDL ON"
)

# seqan3
set (SEQAN3_VERSION 3.4.0 CACHE STRING "")
CPMDeclarePackage (seqan3
                   NAME seqan3
                   VERSION ${SEQAN3_VERSION}
                   GIT_TAG ${SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# fmt
set (FMT_VERSION e424e3f2e607da02742f73db84873b8084fc714c CACHE STRING "")
CPMDeclarePackage (fmt
                   NAME fmt
                   GIT_TAG ${FMT_VERSION} # master
                   GITHUB_REPOSITORY fmtlib/fmt
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING" "FMT_INSTALL OFF"
)

# fmindex
set (FMINDEX_VERSION e35626ba983dde88107709345f1efe49f1ac02e1 CACHE STRING "")
CPMDeclarePackage (fmindex
                   NAME fmindex
                   GIT_TAG ${FMINDEX_VERSION} # main
                   GITHUB_REPOSITORY SGSSGene/fmindex-collection
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# googletest
set (GOOGLETEST_VERSION 1.17.0 CACHE STRING "")
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37 CACHE STRING "")
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)

# thresholding
CPMDeclarePackage (thresholding
                   NAME thresholding
                   URL "${PROJECT_SOURCE_DIR}/contrib/threshold.tar.gz" # main
                   URL_HASH SHA256=4990c7fb9778a2fb8a19794b966d57496ca77bcd708b4cee3c93eea6e5b67d80
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)

# cmake-format: on
