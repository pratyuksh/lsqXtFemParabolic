#include <gtest/gtest.h>

#include "mfem.hpp"
#include <iostream>

#include "../src/sparse_heat/solver.hpp"
#include "../src/sparse_heat/observer.hpp"

using namespace mfem;


std::tuple<int, double, double>
runSolverAndVisualize(sparseHeat::Solver& solver,
                      sparseHeat::Observer& observer)
{
    solver.run();
    auto solutionHandler = solver.getSolutionHandler();
    auto testCase = solver.getTestCase();
    auto disc = solver.getDiscretisation();

    observer.set(testCase, disc);
    observer.visualizeSolutionAtEndTime(*solutionHandler);

    double htMax = solver.getMeshwidthOfFinestTemporalMesh();
    double hxMax = solver.getMeshwidthOfFinestSpatialMesh();

    return {solver.getNumDofs(), htMax, hxMax};
}


/**
 * @brief Tests the least-squares space-time H1-Hdiv discretisation with unitSquare_test1
 */
TEST(SparseSolver, buildAndRunUnitSquareTest1)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/"
              "sparse_heat_solver/sparseHeat_unitSquare_test1.json";
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

    std::string meshDir = "../tests/input/sparse_heat_solver";
    bool loadInitMesh = true;

    sparseHeat::Solver solver(config,
                              testCase,
                              meshDir,
                              numLevels,
                              minSpatialLevel,
                              minTemporalLevel,
                              loadInitMesh);

    sparseHeat::Observer observer(config, numLevels, minTemporalLevel);

    int numDofs;
    double htMax, hxMax;
    std::tie(numDofs, htMax, hxMax)
            = runSolverAndVisualize(solver, observer);

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


TEST(SparseSolver, buildAndRunUnitSquareTest3)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/"
              "sparse_heat_solver/sparseHeat_unitSquare_test3.json";
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

    std::string meshDir = "../tests/input/sparse_heat_solver";
    bool loadInitMesh = true;

    sparseHeat::Solver solver(config,
                              testCase,
                              meshDir,
                              numLevels,
                              minSpatialLevel,
                              minTemporalLevel,
                              loadInitMesh);

    sparseHeat::Observer observer(config, numLevels, minTemporalLevel);

    int numDofs;
    double htMax, hxMax;
    std::tie(numDofs, htMax, hxMax)
            = runSolverAndVisualize(solver, observer);

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
/*
TEST(SparseSolver, buildAndRunLShapedTest1)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/"
              "sparse_heat_solver/sparseHeat_lShaped_test1.json";
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

    std::string meshDir = "../meshes/lShaped/br/deg1";
    bool loadInitMesh = false;

    sparseHeat::Solver solver(config,
                              testCase,
                              meshDir,
                              numLevels,
                              minSpatialLevel,
                              minTemporalLevel,
                              loadInitMesh);

    sparseHeat::Observer observer(config, numLevels, minTemporalLevel);

    int numDofs;
    double htMax, hxMax;
    std::tie(numDofs, htMax, hxMax)
            = runSolverAndVisualize(solver, observer);

    std::cout << htMax << "\t"
              << hxMax << "\t"
              << numDofs << std::endl;
    solver.getDiscretisation()
            ->getSpatialNestedFEHierarchyForTemperature()
            ->getNumDims().Print();
    solver.getDiscretisation()
            ->getSpatialNestedFEHierarchyForHeatFlux()
            ->getNumDims().Print();
}*/
