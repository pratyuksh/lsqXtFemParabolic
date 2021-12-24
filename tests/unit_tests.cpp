#include <gtest/gtest.h>


//! Calls all the unit tests
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

// End of file
