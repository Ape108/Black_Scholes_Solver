/*
Given a Tridagonal Matrix A of shape n X n, where: 

                A = [[a_1, b_1,                  ],
                     [c_2, a_2, b_2,             ],
                     [   , c_3, a_3, b_3,        ],
                        ...
                     [        c_n-1, a_n-1, b_n-1],
                     [                   c_n, a_n]]

We use LU Decomposition to separate the tridiagonal matrix
into a lower diagonal matrix and an upper diagonal matrix: 

                A = LU

        where:

                L = [[1,               ],
                     [l_2, 1,          ],
                     [    l_3, 1,      ],
                        ...
                     [           l_n, 1]],
                
        and 
                U = [[u_1, b_1,        ],
                     [    u_2, b_2,    ],
                        ...
                     [     u_n-1, b_n-1],
                     [              u_n]]


Given: b = [b_1, b_2, ..., b_n]^T,

Since we know L is a lower-diagonal matrix with 1's on the main diagonal,
we can just store all of the information in a vector l.

We will evaluate the linear system: Ly = b
because we will use the values of y to solve Ux = y.

Forward substitution on an upper diagonal matrix of this
specific form is an O(n) operation.

Given: y = [y_1, y_2, ... y_n]^T,

We know that the non-zero values of the upper-diagonal matrix U
lie on the main and upper diagonals, so we can flatten them into 
two vectors containing all of the information we need.

We obtained the y vector from solving the linear system Ly = b,
and now we will use it to solve Ux = y.

x will be the solution vector to the original tridiagonal linear system.

*/

#pragma once

#include <vector>
#include <stdexcept>
#include <string>
#include <cmath>

// Struct to hold the values of L and U. 
struct Decomposed {
    std::vector<double> lower;
    std::vector<double> upper;
};

/// @brief Performs LU decomposition on a tridiagonal matrix.
/// @param a The input matrix main-diagonal values (Flattened into a 1D vector).
/// @param b The input matrix upper-diagonal values (Flattened into a 1D vector).
/// @param c The input matrix lower-diagonal values (Flattened into a 1D vector).
/// @return A 'Decomposed' struct containing an upper and lower matrix.
Decomposed lu_decomposition(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);

/// @brief Solves a system of linear equations using forward substitution.
/// @param lower The lower triangular matrix (flattened into a 1D vector).
/// @param b The right-hand side constant vector.
/// @param y A pre-allocated buffer where the resulting solution vector will be written.
void forward_substitution(const std::vector<double>& lower, const std::vector<double>& b, std::vector<double>& y);

/// @brief Solves a system of linear equations using backward substitution.
/// @param u The main diagonal of the upper triangular matrix (flattened into a 1D vector).
/// @param b The upper diagonal of the upper triangular matrix (flattened into a 1D vector).
/// @param y The right-hand side constant vector.
/// @param x A pre-allocated buffer where the final solution vector will be written.
void backward_substitution(const std::vector<double>& u, const std::vector<double>& b, const std::vector<double>& y, std::vector<double>& x);