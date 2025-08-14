// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <mutex>
#include <string>
#include <vector>

#include <fpgalign/config.hpp>
#include <fpgalign/meta.hpp>

namespace search
{

struct hit
{
    std::vector<uint64_t> bins;
};

struct wip_alignment
{
    size_t bin;
    size_t sequence_number;
    size_t position;
    size_t idx;
};

class alignment_vector
{
private:
    std::vector<wip_alignment> data;
    std::mutex mtx;

public:
    std::vector<wip_alignment> & get() noexcept
    {
        return data;
    }

    std::vector<wip_alignment> const & get() const noexcept
    {
        return data;
    }

    void emplace_back(wip_alignment elem)
    {
        std::lock_guard guard{mtx};
        data.emplace_back(std::move(elem));
    }

    void emplace_back(wip_alignment & elem) = delete;
};

void search(config const & config);
std::vector<hit> ibf(config const & config, meta & meta);
std::vector<wip_alignment> fmindex(config const & config, meta & meta, std::vector<hit> hits);
void do_alignment(config const & config, meta & meta, std::vector<wip_alignment> const & wips);

} // namespace search
