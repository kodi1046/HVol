
#include <numbers>
#include <cmath>
#include <stdexcept>
#include <ranges>


constexpr double pi = std::numbers::pi;
namespace fft {
    namespace detail {
        template <std::ranges::random_access_range R>
        constexpr auto size(R&& r) noexcept {
            return std::ranges::size(std::forward<R>(r));
        }
        
        // radix_2_DIT implementation (default)
        //
        // The Cooley-Tukey radix_2_DIT algorithm is in the family of FFT algorithms. 
        // It separately computes the DFTs of the even and odd indexed inputs, and then combines them.
        // This can be implemented recursively, and brings the runtime down to O(nlogn). 
        // 
        template <std::ranges::random_access_range Range>
            requires ComplexLike<std::ranges::range_value_t<Range>>
        auto fft_forward_impl(Range&& r, tag::radix_2_DIT) {
            using value_type = std::ranges::range_value_t<Range>;
            auto result = std::ranges::to<std::remove_cvref_t<Range>>(std::forward<Range>(r));
            
            std::size_t n = size(result);
            if (n <= 1) return result;
            
            if  ((n & (n - 1)) != 0) {
                throw std::invalid_argument("FFT size must be a power of 2");       
            }
            
            // bit reversal
            std::size_t log2n = std::bit_width(n) - 1;
            
            for (std::size_t i = 1; i < n; ++i) {
                std::size_t rev = utils::bit_reverse(i, log2n);
                if (i < rev) {
                    std::ranges::swap(result[i], result[rev]);
                }
            }
            
            for (std::size_t s = 1; s <= log2n; ++s) {
                std::size_t m = size_t{1} << s;
                std::size_t m2 = m >> 1;
                
                value_type wm = std::polar(1.0, -2 * pi / double(m));
                
                for (std::size_t k = 0; k < n; k += m) {
                    value_type w{1.0};
                    
                    for (std::size_t j = 0; j < m2; ++j) {
                        value_type t = w * result[k + j + m2];
                        value_type u = result[k + j];
                        result[k + j] = u + t;
                        result[k + j + m2] = u - t;
                        w *= wm;
                    }
                }
            
            }
            return result;
        }
            
        template <std::ranges::random_access_range Range>
            requires ComplexLike<std::ranges::range_value_t<Range>>
        auto fft_inverse_impl(Range&& r, tag::radix_2_DIT) {
            using value_type = std::ranges::range_value_t<Range>;
            auto result = std::ranges::to<std::remove_cvref_t<Range>>(std::forward<Range>(r));
            
            for (auto& v : result) v = std::conj(v);
            
            result = fft_forward_impl(std::move(result), tag::radix_2_DIT{});
            
            for (auto& v : result) v = std::conj(v);
            
            auto n = static_cast<typename value_type::value_type>(result.size());
            for (auto& v : result) v /= n;
            
            return result;
        }   
        
        // radix_2_DIF implementation
        template <std::ranges::random_access_range Range>
            requires ComplexLike<std::ranges::range_value_t<Range>>
        auto fft_forward_impl(Range&& r, tag::radix_2_DIF) {
            using value_type = std::ranges::range_value_t<Range>;
            auto result = std::ranges::to<std::remove_cvref_t<Range>>(std::forward<Range>(r));
            
            // TODO: implement
            
            return result;
        }
        
        template <std::ranges::random_access_range Range>
            requires ComplexLike<std::ranges::range_value_t<Range>>
        auto fft_inverse_impl(Range&& r, tag::radix_2_DIF) {
            using value_type = std::ranges::range_value_t<Range>;
            auto result = std::ranges::to<std::remove_cvref_t<Range>>(std::forward<Range>(r));
            
            // TODO: implement
            
            return result;
        }  
    }
}