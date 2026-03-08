#pragma once

#include <concepts>
#include <ranges>
#include <utility>
#include <complex> // IWYU pragma: keep

namespace fft {
     
    template <typename T>
    concept ComplexLike = 
        std::regular<T> &&
        requires(T a, T b) {
            { a + b } -> std::same_as<T>;
            { a - b } -> std::same_as<T>;
            { a * b } -> std::same_as<T>;
        };
    
    namespace tag {
        struct radix_2_DIT {};
        struct radix_2_DIF {};
        // ... add more tags as needed
    }
    
    namespace detail {
        
        // radix_2_DIT
        template <std::ranges::random_access_range R>
            requires ComplexLike<std::ranges::range_value_t<R>>
        auto fft_forward_impl(R&& r, tag::radix_2_DIT);
        
        template <std::ranges::random_access_range R>
            requires ComplexLike<std::ranges::range_value_t<R>>
        auto fft_inverse_impl(R&& r, tag::radix_2_DIT);
        
        // radix_2_DIF
        template <std::ranges::random_access_range R>
            requires ComplexLike<std::ranges::range_value_t<R>>
        auto fft_forward_impl(R&& r, tag::radix_2_DIF);
        
        template <std::ranges::random_access_range R>
            requires ComplexLike<std::ranges::range_value_t<R>>
        auto fft_inverse_impl(R&& r, tag::radix_2_DIF);
        
    }
    
    // Generic forward/backward FFT declaration using tag dispatch.
    // With radix_2_DIT as the default.
    template <std::ranges::random_access_range R,
              typename Tag = tag::radix_2_DIT>
        requires ComplexLike<std::ranges::range_value_t<R>>  
    auto fft_forward(R&& r) {
        return detail::fft_forward_impl(std::forward<R>(r), Tag{});
    }
    
    template <std::ranges::random_access_range R,
              typename Tag = tag::radix_2_DIT>
        requires ComplexLike<std::ranges::range_value_t<R>>
    auto fft_inverse(R&& r) {
        return detail::fft_inverse_impl(std::forward<R>(r), Tag{});
    }
    
    namespace radix_2_DIT {
        using fft::fft_forward;
        using fft::fft_inverse;
    }
    
    namespace radix_2_DIF {
        using fft::fft_forward;
        using fft::fft_inverse;
    }
    
    // ... add more as needed
    
}

// Actual implementation
#include "fft.tpp"

