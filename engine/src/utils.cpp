#include "../include/utils.hpp"

std::size_t utils::bit_reverse(std::size_t x, std::size_t log2n) {
    std::size_t rev = 0;
    for (std::size_t i = 0; i < log2n; ++i) {
        if (x & (std::size_t{1} << i)) {
            rev |= std::size_t{1} << (log2n - 1 - i);
        }
    }
    
    return rev;
}
