#include "linear_algebra.hpp"

Decomposed lu_decomposition(const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& c) {
    
    Decomposed LU;
    LU.upper.push_back(a[0]);
    size_t n = a.size();

    for (size_t i=1; i<n; i++) {
        float l_i = c[i-1] / LU.upper[i-1];
        LU.lower.push_back(l_i);
        float u_i = a[i] - (l_i * b[i-1]);
        LU.upper.push_back(u_i); 
    }

    return LU;
}

std::vector<float> forward_substitution(const std::vector<float>& lower, const std::vector<float>& b) {
    
    std::vector<float> y = {b[0]};
    size_t n = lower.size();
    
    for (size_t i=1; i<=n; i++) {
        float y_i = b[i] - (lower[i-1] * y[i-1]);
        y.push_back(y_i);
    }

    return y;
}

std::vector<float> backward_substitution(const std::vector<float>& u, const std::vector<float>& b, const std::vector<float>& y) {

    size_t n = u.size();
    std::stack<float> x;

    x.push(y[n-1] / u[n-1]);
    
    for (size_t i = n-1; i > 0; i--) {
        float x_i = y[i-1] - (b[i-1] * x.top());
        x_i /= u[i-1];
        x.push(x_i);
    }

    std::vector<float> solution;

    // Pop all values from stack into the vector in correct order
    while (!x.empty()) {
        float top = x.top();
        solution.push_back(top);
        x.pop();
    }

    return solution;
}