#include "../src/main.cpp"

#include "gtest/gtest.h"

TEST(Shadows, triangle1)
{
    EXPECT_EQ(2, 1);
    // EXPECT_EQ(hill_shading(1,0,0,0,1,0,0,0,1,45,45), .5);
}

int sideEffect = 42;
bool f() {
    sideEffect = 16;
    return false;
}

TEST(ExampleTests, DemonstrateGTestMacros)
{
    EXPECT_TRUE(false);
    ASSERT_TRUE(false);
    const bool result = f();
    EXPECT_EQ(true, result) << "sideEffect = " << sideEffect;
}