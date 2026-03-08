#include <numbers>

namespace fft::detail {
    constexpr auto pi = std::numbers::pi;
    
    template <std::ranges::range R>
    constexpr auto size(R&& r) noexcept {
        return std::ranges::size(std::forward<R>(r));
    }
    
    // radix_2_DIT implementation (default)
    //
    // The Cooley-Tukey radix_2_DIT algorithm is in the family of FFT algorithms. 
    // It separately computes the DFTs of the even and odd indexed inputs, and then combines them.
    // This can be implemented recursively, and brings the runtime down to O(nlogn). 
    // 
    template <std::ranges::random_access_range R>
        requires ComplexLike<std::ranges::range_value_t<R>>
    auto fft_forward_impl(R&& r, tag::radix_2_DIT) {
        using value_type = std::ranges::range_value_t<R>;
        auto result = std::ranges::to<std::remove_cvref_t<R>>(std::forward<R>(r));
        
        std::size_t n = size(result);
        if (n <= 1) return result;
        
        auto even = fft::fft_forward_impl(result | std::views::filter([](auto i) { return i % 2 == 0; }), tag::radix_2_DIT);
        auto odd = fft::fft_forward_impl(result | std::views::filter([](auto i) { return i % 2 == 1; }), tag::radix_2_DIT);
        
        for (std::size_t k = 0; k < n; ++k) {
            auto t = std::polar(1.0, -2 * pi * k / n) * odd[k];
            result[k] = even[k] + t;
            result[k + n / 2] = even[k] - t;
        }
        
        return result;
    }
    
    template <std::ranges::random_access_range R>
        requires ComplexLike<std::ranges::range_value_t<R>>
    auto fft_inverse_impl(R&& r, tag::radix_2_DIT) {
        using value_type = std::ranges::range_value_t<R>;
        auto result = std::ranges::to<std::remove_cvref_t<R>>(std::forward<R>(r));
        
        result = std::conj(result); 
        
        fft_forward_impl(result, tag::radix_2_DIT);
        
        result = std::conj(result);
        
        result /= size(result);
        
        return result;
    }
    
    // radix_2_DIF implementation
    template <std::ranges::random_access_range R>
        requires ComplexLike<std::ranges::range_value_t<R>>
    auto fft_forward_impl(R&& r, tag::radix_2_DIF) {
        using value_type = std::ranges::range_value_t<R>;
        auto result = std::ranges::to<std::remove_cvref_t<R>>(std::forward<R>(r));
        
        // TODO: implement
        
        return result;
    }
    
    template <std::ranges::random_access_range R>
        requires ComplexLike<std::ranges::range_value_t<R>>
    auto fft_inverse_impl(R&& r, tag::radix_2_DIF) {
        using value_type = std::ranges::range_value_t<R>;
        auto result = std::ranges::to<std::remove_cvref_t<R>>(std::forward<R>(r));
        
        // TODO: implement
        
        return result;
    }  
}