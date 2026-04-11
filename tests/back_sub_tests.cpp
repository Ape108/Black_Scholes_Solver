#include "linear_algebra.hpp"

#include <iostream>
#include <gtest/gtest.h>

TEST(BackSubTest, BasicAssertions) {
    std::vector<float> u = {1.0,5.0,0.2};
    std::vector<float> b = {1.0,8.0};
    std::vector<float> y = {6.0, -3.0, 7.8};
    std::vector<float> x = backward_substitution(u, b, y);


    EXPECT_NEAR(x[0], 69.0f, 1e-6);
    EXPECT_NEAR(x[1],-63.0f, 1e-6);
    EXPECT_NEAR(x[2], 39.0f, 1e-6);
}

// Test case with single equation (size 1)
TEST(BackSubTest, SingleElement) {
    std::vector<float> u = {2.0};
    std::vector<float> b = {};
    std::vector<float> y = {10.0};
    std::vector<float> x = backward_substitution(u, b, y);

    EXPECT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 5.0f, 1e-6);
}

// Test case with two equations
TEST(BackSubTest, TwoElements) {
    std::vector<float> u = {2.0, 3.0};
    std::vector<float> b = {4.0};
    std::vector<float> y = {12.0, 9.0};
    std::vector<float> x = backward_substitution(u, b, y);

    EXPECT_EQ(x.size(), 2);
    EXPECT_NEAR(x[0], 0.0f, 1e-6);
    EXPECT_NEAR(x[1], 3.0f, 1e-6);
}

// Test with negative values
TEST(BackSubTest, NegativeValues) {
    std::vector<float> u = {-1.0, 2.0, -0.5};
    std::vector<float> b = {3.0, -2.0};
    std::vector<float> y = {-5.0, 4.0, 1.0};
    std::vector<float> x = backward_substitution(u, b, y);

    EXPECT_EQ(x.size(), 3);
    // Verify by checking system is solved correctly
    float check_x2 = y[2] / u[2];
    float check_x1 = (y[1] - b[1] * check_x2) / u[1];
    EXPECT_NEAR(x[0], (y[0] - b[0] * check_x1) / u[0], 1e-5);
    EXPECT_NEAR(x[1], check_x1, 1e-5);
    EXPECT_NEAR(x[2], check_x2, 1e-5);
}

// Test with larger system (size 4)
TEST(BackSubTest, FourElements) {
    std::vector<float> u = {1.0, 2.0, 3.0, 4.0};
    std::vector<float> b = {1.0, 1.0, 1.0};
    std::vector<float> y = {10.0, 20.0, 30.0, 40.0};
    std::vector<float> x = backward_substitution(u, b, y);

    EXPECT_EQ(x.size(), 4);
    float check_x3 = y[3] / u[3];
    float check_x2 = (y[2] - b[2] * check_x3) / u[2];
    float check_x1 = (y[1] - b[1] * check_x2) / u[1];
    EXPECT_NEAR(x[0], (y[0] - b[0] * check_x1) / u[0], 1e-5);
    EXPECT_NEAR(x[1], check_x1, 1e-5);
    EXPECT_NEAR(x[2], check_x2, 1e-5);
    EXPECT_NEAR(x[3], check_x3, 1e-5);
}

// Test with zeros in y vector
TEST(BackSubTest, ZerosInY) {
    std::vector<float> u = {1.0, 1.0, 1.0};
    std::vector<float> b = {0.5, 0.5};
    std::vector<float> y = {2.0, 0.0, 0.0};
    std::vector<float> x = backward_substitution(u, b, y);

    EXPECT_EQ(x.size(), 3);
    EXPECT_NEAR(x[2], 0.0f, 1e-6);
}

// Test with fractional coefficients
TEST(BackSubTest, FractionalCoefficients) {
    std::vector<float> u = {0.5, 0.25, 0.125};
    std::vector<float> b = {0.5, 0.25};
    std::vector<float> y = {1.0, 0.5, 0.25};
    std::vector<float> x = backward_substitution(u, b, y);

    EXPECT_EQ(x.size(), 3);
    float check_x2 = y[2] / u[2];
    float check_x1 = (y[1] - b[1] * check_x2) / u[1];
    EXPECT_NEAR(x[0], (y[0] - b[0] * check_x1) / u[0], 1e-5);
    EXPECT_NEAR(x[1], check_x1, 1e-5);
    EXPECT_NEAR(x[2], check_x2, 1e-5);
}