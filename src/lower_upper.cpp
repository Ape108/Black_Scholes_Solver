#include "lower_upper.hpp"

Decomposed Decomp(const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& c) {
    Decomposed LU;
    LU.upper.push_back(a[0]); // Set u_1 <--- a_1
    std::size_t n = a.size();

    for (std::size_t i=1; i<n; i++) {
        float l_i = c[i-1] / LU.upper[i-1]; // l_i = c_i / u_i-1
        LU.lower.push_back(l_i);
        float u_i = a[i] - (l_i * b[i-1]);
        LU.upper.push_back(u_i); 
    }

    return LU;
}
