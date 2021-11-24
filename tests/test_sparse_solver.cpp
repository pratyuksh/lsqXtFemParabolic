#include <gtest/gtest.h>

#include "mfem.hpp"
#include <iostream>

#include "../src/sparse_heat/solver.hpp"

using namespace mfem;

/**
 * @brief Tests the least-squares space-time H1-Hdiv discretisation
 */
TEST(SparseSolver, buildComponentsAndRun)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/sparse_solver/sparseHeat_square_test1.json";
    auto config = getGlobalConfig(configFile);
    auto testCase = heat::makeTestCase(config);

    int deg = config["deg"];
    ASSERT_EQ(deg, 1);

    int numLevels = 1;
    if (config.contains("num_levels")) {
        numLevels = config["num_levels"];
    }
    int minSpatialLevel = 1;
    if (config.contains("min_spatial_level")) {
        minSpatialLevel = config["min_spatial_level"];
    }

    int minTemporalLevel = 1;
    if (config.contains("min_temporal_level")) {
        minTemporalLevel = config["min_temporal_level"];
    }

    std::string meshDir = "../tests/input/sparse_solver";
    bool loadInitMesh = true;

    sparseHeat::Solver solver(config,
                              testCase,
                              meshDir,
                              numLevels,
                              minSpatialLevel,
                              minTemporalLevel,
                              loadInitMesh);

    int numDofs;
    double htMax, hxMax;
    std::tie(numDofs, htMax, hxMax, std::ignore) = solver();

    std::cout << htMax << "\t"
              << hxMax << "\t"
              << numDofs << std::endl;
    solver.getDiscretisation()
            ->getSpatialNestedFEHierarchyForTemperature()
            ->getNumDims().Print();
    solver.getDiscretisation()
            ->getSpatialNestedFEHierarchyForHeatFlux()
            ->getNumDims().Print();
}
