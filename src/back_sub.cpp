#include "../include/back_sub.hpp"

#include <stack>

std::vector<float> backward_substitution(const std::vector<float>& u, const std::vector<float>& b, const std::vector<float>& y) {

    size_t n = u.size();
    std::stack<float> x;

    x.push(y[n] / u[n]);
    
    // i-- > 0 avoids size_t becoming negative, ex: (size_t i=n-1; i >= 0; i--)
    for (size_t i = n; i-- > 0;) {
        float x_i = y[i] - (b[i+1] * y[i+1]);
        x_i = x_i / u[i];
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