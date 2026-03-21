#include "../include/lower_upper.hpp"
#include "../include/forward_sub.hpp"
#include "../include/back_sub.hpp"

#include <iostream>
#include <iomanip>

struct TridiagonalMatrix {
    std::vector<float> main_diag;
    std::vector<float> upper_diag;
    std::vector<float> lower_diag;
};

/// @brief Algorithm for solving tridiagonal linear systems with O(n) time complexity.
/// @param A Tridiagonal matrix.
/// @param rhs Right-hand-side vector.
/// @return Solution vector to the linear system.
std::vector<float> thomas_algorithm(const TridiagonalMatrix& A, const std::vector<float>& rhs) {

    Decomposed LU_matrices = lu_decomposition(A.main_diag, A.upper_diag, A.lower_diag);

    std::vector<float> lower_matrix = LU_matrices.lower;
    std::vector<float> y = forward_substitution(lower_matrix, rhs); // Solution Vector Ly = rhs

    std::vector<float> upper_matrix = LU_matrices.upper;
    std::vector<float> x = backward_substitution(upper_matrix, A.upper_diag, y); // Solution Vector Ux = y

    return x;
}

template <typename T>
void print_vector(const std::vector<T>& vec) {
    size_t n = vec.size();
    for (size_t i=0; i<n; i++) {
        std::cout << std::fixed << std::setprecision(2) << vec[i] << std::endl;
    }
    
}

int main() {

    // A = [1 1 0]
    //     [2 7 8]
    //     [0 3 5]
    TridiagonalMatrix A{{1, 7, 5}, {1, 8}, {2, 3}};

    // Right-hand side
    std::vector<float> rhs = {6.0, 9.0, 6.0};

    // Solution Vector
    std::vector<float> x = thomas_algorithm(A, rhs);

    print_vector(x);

    return 0;
}