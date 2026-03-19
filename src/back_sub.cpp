#include "../include/back_sub.hpp"

#include <stack>

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