#include "../include/forward_sub.hpp"

std::vector<float> forward_substitution(const std::vector<float>& lower, const std::vector<float>& b) {
    
    std::vector<float> y = {b[0]};
    size_t n = lower.size();
    
    for (size_t i=1; i<=n; i++) {
        float y_i = b[i] - (lower[i-1] * y[i-1]);
        y.push_back(y_i);
    }

    return y;
}


