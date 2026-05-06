#include "linear_algebra.hpp"

Decomposed lu_decomposition(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c) {
    
    if (a.empty()) {
        throw std::invalid_argument("LU decomposition: main diagonal vector 'a' cannot be empty");
    }
    
    size_t n = a.size();
    
    if (std::abs(a[0]) < 1e-15) {
        throw std::runtime_error("LU decomposition: first pivot element is zero or near-zero, matrix is singular");
    }
    
    Decomposed LU;
    LU.upper.push_back(a[0]);

    for (size_t i=1; i<n; i++) {
        double denominator = LU.upper[i-1];
        if (std::abs(denominator) < 1e-15) {
            throw std::runtime_error("LU decomposition: pivot became zero or near-zero at iteration " + std::to_string(i) + ", matrix is singular");
        }
        
        double l_i = c[i-1] / denominator;
        LU.lower.push_back(l_i);
        double u_i = a[i] - (l_i * b[i-1]);
        LU.upper.push_back(u_i); 
    }

    return LU;
}

void forward_substitution(const std::vector<double>& lower, const std::vector<double>& b, std::vector<double>& y) {
    
    size_t n = lower.size();
    y[0] = b[0];

    for (size_t i=1; i<=n; i++) {
        y[i] = b[i] - (lower[i-1] * y[i-1]);
    }
}

void backward_substitution(const std::vector<double>& u, const std::vector<double>& b, const std::vector<double>& y, std::vector<double>& x) {
    
    size_t n = u.size();

    x[n - 1] = y[n - 1] / u[n - 1];

    for (size_t i = n - 1; i > 0; --i) {
        x[i - 1] = (y[i - 1] - (b[i - 1] * x[i])) / u[i - 1];
    }
}