#include "../include/forward_sub.hpp"

#include <iostream>
#include <gtest/gtest.h>

TEST(ForwardSubTest, BasicAssertions) {
    std::vector<float> l = {2.0, 0.6};
    std::vector<float> b = {6.0, 9.0, 6.0};
    std::vector<float> y = forward_substitution(l, b);


    EXPECT_NEAR(y[0], 6.0f, 1e-6);
    EXPECT_NEAR(y[1],-3.0f, 1e-6);
    EXPECT_NEAR(y[2], 7.8f, 1e-6);
}

// Test case with single element
TEST(ForwardSubTest, SingleElement) {
    std::vector<float> l = {};
    std::vector<float> b = {5.0};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 1);
    EXPECT_NEAR(y[0], 5.0f, 1e-6);
}

// Test case with two elements
TEST(ForwardSubTest, TwoElements) {
    std::vector<float> l = {2.0};
    std::vector<float> b = {4.0, 8.0};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 2);
    EXPECT_NEAR(y[0], 4.0f, 1e-6);
    EXPECT_NEAR(y[1], 0.0f, 1e-6);
}

// Test with negative values
TEST(ForwardSubTest, NegativeValues) {
    std::vector<float> l = {-1.0, 0.5};
    std::vector<float> b = {10.0, -5.0, 8.0};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 10.0f, 1e-6);
    EXPECT_NEAR(y[1], (-5.0 - (-1.0) * 10.0), 1e-6);
    float expected_y2 = 8.0 - 0.5 * y[1];
    EXPECT_NEAR(y[2], expected_y2, 1e-5);
}

// Test with larger system (size 4)
TEST(ForwardSubTest, FourElements) {
    std::vector<float> l = {1.0, 2.0, 3.0};
    std::vector<float> b = {5.0, 10.0, 15.0, 20.0};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 4);
    EXPECT_NEAR(y[0], 5.0f, 1e-6);
    EXPECT_NEAR(y[1], b[1] - l[0] * y[0], 1e-6);
    EXPECT_NEAR(y[2], b[2] - l[1] * y[1], 1e-6);
    EXPECT_NEAR(y[3], b[3] - l[2] * y[2], 1e-6);
}

// Test with zeros in lower vector
TEST(ForwardSubTest, ZerosInLower) {
    std::vector<float> l = {0.0, 0.0};
    std::vector<float> b = {3.0, 5.0, 7.0};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 3.0f, 1e-6);
    EXPECT_NEAR(y[1], 5.0f, 1e-6);
    EXPECT_NEAR(y[2], 7.0f, 1e-6);
}

// Test with fractional values
TEST(ForwardSubTest, FractionalValues) {
    std::vector<float> l = {0.5, 0.25};
    std::vector<float> b = {2.0, 1.0, 0.5};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 2.0f, 1e-6);
    EXPECT_NEAR(y[1], 1.0 - 0.5 * y[0], 1e-6);
    EXPECT_NEAR(y[2], 0.5 - 0.25 * y[1], 1e-6);
}

// Test with larger coefficients
TEST(ForwardSubTest, LargeCoefficients) {
    std::vector<float> l = {10.0, 5.0};
    std::vector<float> b = {100.0, 200.0, 150.0};
    std::vector<float> y = forward_substitution(l, b);

    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 100.0f, 1e-5);
    EXPECT_NEAR(y[1], 200.0 - 10.0 * y[0], 1e-5);
    EXPECT_NEAR(y[2], 150.0 - 5.0 * y[1], 1e-5);
}