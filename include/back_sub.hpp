/*
Given  U = [[u_1, b_1,        ],
            [    u_2, b_2,    ],
            ...
            [     u_n-1, b_n-1],
            [              u_n]],

and    y = [y_1, y_2, ... y_n]^T,

We know that the non-zero values of the upper-diagonal matrix U
lie on the main and upper diagonals, so we can flatten them into 
two vectors containing all of the information we need.

We obtained the y vector from solving the linear system Ly = b,
and now we will use it to solve Ux = y.

x will be the solution vector to the original tridiagonal linear system.

*/


#pragma once

#include <vector>

/// @brief Solves a system of linear equations using backward substitution.
/// @param u The main diagonal of the upper triangular matrix (flattened into a 1D vector).
/// @param b The upper diagonal of the upper triangular matrix (flattened into a 1D vector).
/// @param y The right-hand side constant vector.
/// @return A vector containing the solution to the system.
std::vector<float> backward_substitution(const std::vector<float>& u, const std::vector<float>& b, const std::vector<float>& y);