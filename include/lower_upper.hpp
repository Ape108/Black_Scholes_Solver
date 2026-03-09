#pragma once

#include <vector>

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


*/

