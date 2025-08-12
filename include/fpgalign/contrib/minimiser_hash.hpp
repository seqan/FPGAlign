// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <hibf/contrib/std/detail/adaptor_from_functor.hpp>

// Same as seqan3::views::minimiser_hash, but optimized for dna4, with integrated adjust_seed

namespace contrib
{

struct minimiser_hash_parameters
{
    uint8_t kmer_size{};
    uint32_t window_size{};
};

} // namespace contrib

namespace contrib::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class minimiser_hash : public std::ranges::view_interface<minimiser_hash<range_t>>
{
private:
    range_t range{};
    minimiser_hash_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    minimiser_hash()
        requires std::default_initializable<range_t>
    = default;
    minimiser_hash(minimiser_hash const & rhs) = default;
    minimiser_hash(minimiser_hash && rhs) = default;
    minimiser_hash & operator=(minimiser_hash const & rhs) = default;
    minimiser_hash & operator=(minimiser_hash && rhs) = default;
    ~minimiser_hash() = default;

    explicit minimiser_hash(range_t range, minimiser_hash_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class minimiser_hash<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = uint64_t;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{};
    uint64_t seed{};

    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};
    size_t minimiser_position{};

    size_t range_size{};
    size_t range_position{};

    value_type minimiser_value{};

    int kmer_rev_shift{};

    std::deque<uint64_t> kmer_values_in_window{};

    static inline constexpr uint64_t compute_mask(uint8_t const kmer_size)
    {
        assert(kmer_size > 0u);
        assert(kmer_size <= 32u);

        if (kmer_size == 32u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (2u * kmer_size)) - 1u;
    }

    static inline constexpr uint64_t compute_seed(uint8_t const kmer_size)
    {
        assert(kmer_size > 0u);
        assert(kmer_size <= 32u);

        return uint64_t{0x8F3F73B5CF1C9ADEULL} >> (64u - 2u * kmer_size);
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        seed{it.seed},
        kmer_value{it.kmer_value},
        minimiser_position{it.minimiser_position},
        range_size{it.range_size},
        range_position{it.range_position},
        minimiser_value{it.minimiser_value},
        kmer_rev_shift{it.kmer_rev_shift},
        kmer_values_in_window{it.kmer_values_in_window}
    {}

    basic_iterator(range_iterator_t range_iterator, size_t const range_size, minimiser_hash_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(params.kmer_size)},
        seed{compute_seed(params.kmer_size)},
        range_size{range_size},
        kmer_rev_shift{2 * static_cast<int>(params.kmer_size - 1)}
    {
        if (range_size < params.window_size)
            range_position = range_size;
        else
            init(params);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position == lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser_is_new())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser_is_new())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return minimiser_value;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
        requires std::same_as<std::ranges::range_value_t<range_t>, uint8_t>
    {
        uint64_t const new_rank = *range_it;

        kmer_value <<= 2;
        kmer_value |= new_rank;
        kmer_value &= kmer_mask;

        kmer_value_rev >>= 2;
        kmer_value_rev |= (new_rank ^ 0b11) << kmer_rev_shift;
    }

    void rolling_hash()
        requires std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>
    {
        uint64_t const new_rank = range_it->to_rank();

        kmer_value <<= 2;
        kmer_value |= new_rank;
        kmer_value &= kmer_mask;

        kmer_value_rev >>= 2;
        kmer_value_rev |= (new_rank ^ 0b11) << kmer_rev_shift;
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes)
            kmer_values_in_window.pop_front();

        kmer_values_in_window.push_back(std::min<uint64_t>(kmer_value ^ seed, kmer_value_rev ^ seed));
    }

    void find_minimiser_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values_in_window, std::less_equal<uint64_t>{});
        minimiser_value = *minimiser_it;
        minimiser_position = std::distance(std::begin(kmer_values_in_window), minimiser_it);
    }

    void init(minimiser_hash_parameters const & params)
    {
        // range_it is already at the beginning of the range
        rolling_hash();
        kmer_values_in_window.push_back(std::min<uint64_t>(kmer_value ^ seed, kmer_value_rev ^ seed));

        // After this loop, `kmer_values_in_window` contains the first kmer value of the window.
        for (size_t i = 1u; i < params.kmer_size; ++i)
            next_window<pop_first::yes>();

        // After this loop, `kmer_values_in_window` contains all kmer values of the window.
        for (size_t i = params.kmer_size; i < params.window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser_in_window();
    }

    bool next_minimiser_is_new()
    {
        // If we reached the end of the range, we are done.
        if (range_position + 1 == range_size)
            return ++range_position; // Return true, but also increment range_position

        next_window<pop_first::yes>();

        // The minimiser left the window.
        if (minimiser_position == 0)
        {
            find_minimiser_in_window();
            return true;
        }

        // Update minimiser if the new kmer value is smaller than the current minimiser.
        if (uint64_t new_kmer_value = kmer_values_in_window.back(); new_kmer_value < minimiser_value)
        {
            minimiser_value = new_kmer_value;
            minimiser_position = kmer_values_in_window.size() - 1u;
            return true;
        }

        --minimiser_position;
        return false;
    }
};

template <std::ranges::viewable_range rng_t>
minimiser_hash(rng_t &&, minimiser_hash_parameters &&) -> minimiser_hash<std::views::all_t<rng_t>>;

struct minimiser_hash_fn
{
    constexpr auto operator()(minimiser_hash_parameters params) const
    {
        return seqan::stl::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, minimiser_hash_parameters params) const
    {
        // static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.kmer_size == 0u)
            throw std::invalid_argument{"kmer_size must be > 0."};
        if (params.kmer_size > 32u)
            throw std::invalid_argument{"kmer_size must be <= 32."};
        if (params.window_size == 0u)
            throw std::invalid_argument{"window_size must be > 0."};
        if (params.window_size < params.kmer_size)
            throw std::invalid_argument{"window_size must be >= kmer_size."};

        return minimiser_hash{std::forward<range_t>(range), std::move(params)};
    }
};

} // namespace contrib::detail

namespace contrib::views
{

inline constexpr auto minimiser_hash = contrib::detail::minimiser_hash_fn{};

}
