#include "lower_upper.hpp"
#include <iostream>
#include <string>
#include <iterator>

int main() {
    std::vector<float> a = {1, 7, 5};
    std::vector<float> b = {1, 8};
    std::vector<float> c = {2, 3};

    Decomposed LU = Decomp(a, b, c);

    std::cout << "Expected: {2.0, 0.6}" << std::endl;
    std::copy(LU.lower.begin(), LU.lower.end(), std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Expected: {1.0, 5.0, 0.2}" << std::endl;
    std::copy(LU.upper.begin(), LU.upper.end(), std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;
    return 0;
}
