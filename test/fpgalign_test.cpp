// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct fpgalign : public app_test
{};

TEST_F(fpgalign, no_options)
{}
