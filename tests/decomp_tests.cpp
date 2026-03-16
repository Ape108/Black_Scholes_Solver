#include "../include/lower_upper.hpp"

#include <iostream>
#include <gtest/gtest.h>

TEST(DecompTestThree, BasicAssertions) {
    std::vector<float> a = {1, 7, 5};
    std::vector<float> b = {1, 8};
    std::vector<float> c = {2, 3};

    Decomposed LU = Decomp(a, b, c);

    EXPECT_NEAR(LU.lower[0], 2.0f, 1e-6);
    EXPECT_NEAR(LU.lower[1], 0.6f, 1e-6);
    EXPECT_NEAR(LU.upper[0], 1.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1], 5.0f, 1e-6);
    EXPECT_NEAR(LU.upper[2], 0.2f, 1e-6);
}

TEST(DecompTestFour, BasicAssertions) {
    std::vector<float> a = {3, 3, 6, 4};
    std::vector<float> b = {-2, 1, -1};
    std::vector<float> c = {-3, 2, -8};

    Decomposed LU = Decomp(a, b, c);

    EXPECT_NEAR(LU.lower[0],-1.0f, 1e-6);
    EXPECT_NEAR(LU.lower[1], 2.0f, 1e-6);
    EXPECT_NEAR(LU.lower[2],-2.0f, 1e-6);
    EXPECT_NEAR(LU.upper[0], 3.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1], 1.0f, 1e-6);
    EXPECT_NEAR(LU.upper[2], 4.0f, 1e-6);
    EXPECT_NEAR(LU.upper[3], 2.0f, 1e-6);

}

TEST(DecompTestFractions, BasicAssertions) {
    std::vector<float> a = {2, 2, 2, 2};
    std::vector<float> b = {-1, -1, -1};
    std::vector<float> c = {-1, -1, -1};

    Decomposed LU = Decomp(a, b, c);

    EXPECT_NEAR(LU.lower[0],(-1.0/2.0), 1e-6);
    EXPECT_NEAR(LU.lower[1],(-2.0/3.0), 1e-6);
    EXPECT_NEAR(LU.lower[2],(-3.0/4.0), 1e-6);
    EXPECT_NEAR(LU.upper[0], 2.0f, 1e-6);
    EXPECT_NEAR(LU.upper[1],(3.0/2.0), 1e-6);
    EXPECT_NEAR(LU.upper[2],(4.0/3.0), 1e-6);
    EXPECT_NEAR(LU.upper[3],(5.0/4.0), 1e-6);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
