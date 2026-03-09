#include "lower_upper.hpp"
#include <iostream>
#include <string>
#include <iterator>

int main() {
    std::vector<float> a = {1, 7, 5};
    std::vector<float> b = {1, 8};
    std::vector<float> c = {2, 3};

    Decomposed LU = Decomp(a, b, c);
    std::copy(LU.lower.begin(), LU.lower.end(), std::ostream_iterator<float>(std::cout, " "));
    std::copy(LU.upper.begin(), LU.upper.end(), std::ostream_iterator<float>(std::cout, " "));
};
