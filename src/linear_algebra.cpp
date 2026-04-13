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

std::vector<double> forward_substitution(const std::vector<double>& lower, const std::vector<double>& b) {
    
    if (b.empty()) {
        throw std::invalid_argument("Forward substitution: RHS vector 'b' cannot be empty");
    }
    
    std::vector<double> y;
    y.push_back(b[0]);
    size_t n = lower.size();
    
    for (size_t i=1; i<=n; i++) {
        double y_i = b[i] - (lower[i-1] * y[i-1]);
        y.push_back(y_i);
    }

    return y;
}

std::vector<double> backward_substitution(const std::vector<double>& u, const std::vector<double>& b, const std::vector<double>& y) {

    if (u.empty() || y.empty()) {
        throw std::invalid_argument("Backward substitution: upper diagonal 'u' and RHS 'y' cannot be empty");
    }
    
    size_t n = u.size();
    
    if (std::abs(u[n-1]) < 1e-15) {
        throw std::runtime_error("Backward substitution: pivot element u[" + std::to_string(n-1) + "] is zero or near-zero");
    }
    
    std::stack<double> x;
    x.push(y[n-1] / u[n-1]);
    
    for (size_t i = n-1; i > 0; i--) {
        
        if (std::abs(u[i-1]) < 1e-15) {
            throw std::runtime_error("Backward substitution: pivot element u[" + std::to_string(i-1) + "] is zero or near-zero");
        }
        
        double x_i = y[i-1] - (b[i-1] * x.top());
        x_i /= u[i-1];
        x.push(x_i);
    }

    std::vector<double> solution;

    // Pop all values from stack into the vector in correct order
    while (!x.empty()) {
        double top = x.top();
        solution.push_back(top);
        x.pop();
    }

    return solution;
}