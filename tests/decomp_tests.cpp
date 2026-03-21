#include "../include/lower_upper.hpp"

#include <iostream>
#include <gtest/gtest.h>

TEST(DecompTest, BasicAssertions) {
    std::vector<float> a = {1, 7, 5};
    std::vector<float> b = {1, 8};
    std::vector<float> c = {2, 3};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_NEAR(LU.lower[0], 2.0f, 1e-6);
    EXPECT_NEAR(LU.lower[1], 0.6f, 1e-6);
    EXPECT_NEAR(LU.upper[0], 1.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1], 5.0f, 1e-6);
    EXPECT_NEAR(LU.upper[2], 0.2f, 1e-6);
}

TEST(DecompTest, FourElements) {
    std::vector<float> a = {3, 3, 6, 4};
    std::vector<float> b = {-2, 1, -1};
    std::vector<float> c = {-3, 2, -8};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_NEAR(LU.lower[0],-1.0f, 1e-6);
    EXPECT_NEAR(LU.lower[1], 2.0f, 1e-6);
    EXPECT_NEAR(LU.lower[2],-2.0f, 1e-6);
    EXPECT_NEAR(LU.upper[0], 3.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1], 1.0f, 1e-6);
    EXPECT_NEAR(LU.upper[2], 4.0f, 1e-6);
    EXPECT_NEAR(LU.upper[3], 2.0f, 1e-6);

}

TEST(DecompTest, Fractional) {
    std::vector<float> a = {2, 2, 2, 2};
    std::vector<float> b = {-1, -1, -1};
    std::vector<float> c = {-1, -1, -1};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_NEAR(LU.lower[0],(-1.0/2.0), 1e-6);
    EXPECT_NEAR(LU.lower[1],(-2.0/3.0), 1e-6);
    EXPECT_NEAR(LU.lower[2],(-3.0/4.0), 1e-6);
    EXPECT_NEAR(LU.upper[0], 2.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1],(3.0/2.0), 1e-6);
    EXPECT_NEAR(LU.upper[2],(4.0/3.0), 1e-6);
    EXPECT_NEAR(LU.upper[3],(5.0/4.0), 1e-6);
}

// Test with single element (smallest tridiagonal matrix)
TEST(DecompTest, SingleElement) {
    std::vector<float> a = {5.0};
    std::vector<float> b = {};
    std::vector<float> c = {};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 0);
    EXPECT_EQ(LU.upper.size(), 1);
    EXPECT_NEAR(LU.upper[0], 5.0f, 1e-6);
}

// Test with two elements
TEST(DecompTest, TwoElements) {
    std::vector<float> a = {4.0, 3.0};
    std::vector<float> b = {2.0};
    std::vector<float> c = {1.0};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 1);
    EXPECT_EQ(LU.upper.size(), 2);
    EXPECT_NEAR(LU.lower[0], 0.25f, 1e-6);
    EXPECT_NEAR(LU.upper[0], 4.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1], 2.5f, 1e-6);
}

// Test identity-like matrix
TEST(DecompTest, IdentityLike) {
    std::vector<float> a = {1.0, 1.0, 1.0};
    std::vector<float> b = {0.0, 0.0};
    std::vector<float> c = {0.0, 0.0};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 2);
    EXPECT_EQ(LU.upper.size(), 3);
    EXPECT_NEAR(LU.lower[0], 0.0f, 1e-6);
    EXPECT_NEAR(LU.lower[1], 0.0f, 1e-6);
    EXPECT_NEAR(LU.upper[0], 1.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1], 1.0f, 1e-6);
    EXPECT_NEAR(LU.upper[2], 1.0f, 1e-6);
}

// Test with larger system (5 elements)
TEST(DecompTest, FiveElements) {
    std::vector<float> a = {2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<float> b = {1.0, 1.0, 1.0, 1.0};
    std::vector<float> c = {1.0, 1.0, 1.0, 1.0};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 4);
    EXPECT_EQ(LU.upper.size(), 5);
    
    // Verify first decomposition step
    EXPECT_NEAR(LU.upper[0], 2.0f, 1e-6);
    EXPECT_NEAR(LU.lower[0], 1.0f / 2.0f, 1e-6);
}

// Test with fractional main diagonal
TEST(DecompTest, SmallDiagonalValues) {
    std::vector<float> a = {0.5, 0.5, 0.5};
    std::vector<float> b = {0.1, 0.1};
    std::vector<float> c = {0.1, 0.1};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 2);
    EXPECT_EQ(LU.upper.size(), 3);
    
    // Verify first step calculation
    EXPECT_NEAR(LU.upper[0], 0.5f, 1e-6);
    EXPECT_NEAR(LU.lower[0], 0.1f / 0.5f, 1e-6);
}

// Test with negative values
TEST(DecompTest, NegativeValues) {
    std::vector<float> a = {-4.0, 2.0, 3.0};
    std::vector<float> b = {1.0, -1.0};
    std::vector<float> c = {2.0, -1.0};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 2);
    EXPECT_EQ(LU.upper.size(), 3);
    EXPECT_NEAR(LU.upper[0], -4.0f, 1e-6);
    EXPECT_NEAR(LU.lower[0], 2.0f / (-4.0f), 1e-6);
}

// Test with large coefficient variations
TEST(DecompTest, LargeValues) {
    std::vector<float> a = {100.0, 200.0, 150.0, 250.0};
    std::vector<float> b = {50.0, 75.0, 100.0};
    std::vector<float> c = {25.0, 50.0, 75.0};

    Decomposed LU = lu_decomposition(a, b, c);

    EXPECT_EQ(LU.lower.size(), 3);
    EXPECT_EQ(LU.upper.size(), 4);
    
    // Verify first step
    EXPECT_NEAR(LU.upper[0], 100.0f, 1e-6);
    EXPECT_NEAR(LU.lower[0], 25.0f / 100.0f, 1e-6);
}

