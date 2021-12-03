#include <gtest/gtest.h>

#include "mfem.hpp"
#include <iostream>

#include "../src/heat/test_cases_factory.hpp"
#include "../src/heat/discretisation.hpp"
#include "../src/sparse_heat/discretisation.hpp"
#include "../src/sparse_heat/solution_handler.hpp"
#include "../src/sparse_heat/utilities.hpp"

using namespace mfem;

/**
 * @brief Dry run of the least-squares sparse space-time H1-Hdiv discretisation
 */
TEST(SparseDiscretisation, buildH1Hdiv)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/sparse_heat_discretisation/sparseHeat_dummy.json";
    auto config = getGlobalConfig(configFile);
    auto testCase = heat::makeTestCase(config);

    int deg = config["deg"];
    ASSERT_EQ(deg, 1);

    // spatial mesh hierarchy
    std::string input_dir
            = "../tests/input/sparse_heat_assembly/";
    const std::string spatialMeshFile1 = input_dir+"mesh_lx0";
    const std::string spatialMeshFile2 = input_dir+"mesh_lx1";
    const std::string spatialMeshFile3 = input_dir+"mesh_lx2";
    auto spatialMesh1 = std::make_shared<Mesh>(spatialMeshFile1.c_str());
    auto spatialMesh2 = std::make_shared<Mesh>(spatialMeshFile2.c_str());
    auto spatialMesh3 = std::make_shared<Mesh>(spatialMeshFile3.c_str());
    ASSERT_EQ(spatialMesh1->GetNE(), 1);
    ASSERT_EQ(spatialMesh2->GetNE(), 2);
    ASSERT_EQ(spatialMesh3->GetNE(), 3);

    auto spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();
    spatialMeshHierarchy->addMesh(spatialMesh1);
    spatialMeshHierarchy->addMesh(spatialMesh2);
    spatialMeshHierarchy->addMesh(spatialMesh3);
    spatialMeshHierarchy->finalize();

    // temporal hierarchy info
    int minTemporalLevel = 1;
    if (config.contains("min_temporal_level")) {
        minTemporalLevel = config["min_temporal_level"];
    }
    int numLevels = 1;
    if (config.contains("num_levels")) {
        numLevels = config["num_levels"];
    }
    int maxTemporalLevel = minTemporalLevel + numLevels - 1;
    ASSERT_EQ(minTemporalLevel, 1);
    ASSERT_EQ(numLevels, 3);
    ASSERT_EQ(maxTemporalLevel, 3);

    // discretisation
    std::unique_ptr<sparseHeat::LsqSparseXtFem> xtDisc
            = std::make_unique<sparseHeat::LsqSparseXtFemH1Hdiv>
            (config, testCase, numLevels, minTemporalLevel);
    xtDisc->setNestedFEHierarchyAndSpatialBoundaryDofs
            (spatialMeshHierarchy);

    // system assembly
    xtDisc->assembleSystemSubMatrices();
    xtDisc->buildSystemMatrix();

    // rhs assembly
    auto spatialNestedFEHierarchyForTemperature
            = xtDisc->getSpatialNestedFEHierarchyForTemperature();
    auto spatialNestedFEHierarchyForHeatFlux
            = xtDisc->getSpatialNestedFEHierarchyForHeatFlux();

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (minTemporalLevel, maxTemporalLevel,
             spatialNestedFEHierarchyForTemperature,
             spatialNestedFEHierarchyForHeatFlux);

    auto dataSize = solutionHandler->getDataSize();
    auto dataOffsets = evalBlockOffsets(dataSize);
    auto rhs = std::make_shared<BlockVector>(dataOffsets);
    (*rhs) = 0;
    xtDisc->assembleRhs(rhs);

    // verify sizes
    auto systemBlock11 = xtDisc->getSystemBlock11();
    ASSERT_EQ(systemBlock11->NumRows(), 35);
    ASSERT_EQ(systemBlock11->NumCols(), 35);

    auto systemBlock12 = xtDisc->getSystemBlock12();
    ASSERT_EQ(systemBlock12->NumRows(), 35);
    ASSERT_EQ(systemBlock12->NumCols(), 43);

    auto systemBlock21 = xtDisc->getSystemBlock21();
    ASSERT_EQ(systemBlock21->NumRows(), 43);
    ASSERT_EQ(systemBlock21->NumCols(), 35);

    auto systemBlock22 = xtDisc->getSystemBlock22();
    ASSERT_EQ(systemBlock22->NumRows(), 43);
    ASSERT_EQ(systemBlock22->NumCols(), 43);

    auto systemMatrix = xtDisc->getSystemMatrix();
    ASSERT_EQ(systemMatrix->NumRows(), 78);
    ASSERT_EQ(systemMatrix->NumCols(), 78);
}

/**
 * @brief Tests the least-squares sparse space-time H1-Hdiv discretisation with unitSquare_test1
 *
 * Compares the sparse and full discretisations for a single level, with test case unitSquare_test1
 */
TEST(SparseDiscretisation, singleLevelH1HdivUnitSquareTest1)
{
    // config files for full and sparse versions
    std::string heatConfigFile
            = "../config_files/unit_tests/"
              "sparse_heat_discretisation/heat_square_test1.json";
    std::string sparseHeatConfigFile
            = "../config_files/unit_tests/"
              "sparse_heat_discretisation/sparseHeat_square_test1.json";
    auto heatConfig = getGlobalConfig(heatConfigFile);
    auto sparseHeatConfig = getGlobalConfig(sparseHeatConfigFile);

    auto testCase = heat::makeTestCase(heatConfig);

    int deg = heatConfig["deg"];
    ASSERT_EQ(deg, 1);

    int numLevels = 1;
    int spatialLevel = heatConfig["spatial_level"];
    int temporalLevel = heatConfig["temporal_level"];
    ASSERT_EQ(numLevels, sparseHeatConfig["num_levels"]);
    ASSERT_EQ(spatialLevel, sparseHeatConfig["min_spatial_level"]);
    ASSERT_EQ(temporalLevel, sparseHeatConfig["min_temporal_level"]);

    // spatial mesh hierarchy
    std::string inputDir
            = "../tests/input/sparse_heat_discretisation/";
    const std::string spatialMeshFile = inputDir+"mesh_l0.mesh";
    auto spatialMesh
            = std::make_shared<Mesh>(spatialMeshFile.c_str());
    for (int m=0; m<spatialLevel; m++) {
        spatialMesh->UniformRefinement();
    }
    auto spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();
    spatialMeshHierarchy->addMesh(spatialMesh);
    spatialMeshHierarchy->finalize();

    // temporal mesh
    double endTime = heatConfig["end_time"];
    int Nt = static_cast<int>(std::pow(2, temporalLevel));
    auto temporalMesh = std::make_shared<Mesh>(Nt, endTime);

    // sparse system assembly
    std::unique_ptr<sparseHeat::LsqSparseXtFem> sparseXtDisc
            = std::make_unique<sparseHeat::LsqSparseXtFemH1Hdiv>
            (sparseHeatConfig, testCase,
             numLevels, temporalLevel,
             spatialMeshHierarchy);
    sparseXtDisc->assembleSystemSubMatrices();
    sparseXtDisc->buildSystemMatrix();
    auto systemBlock11 = sparseXtDisc->getSystemBlock11();
    auto systemBlock12 = sparseXtDisc->getSystemBlock12();
    auto systemBlock21 = sparseXtDisc->getSystemBlock21();
    auto systemBlock22 = sparseXtDisc->getSystemBlock22();
    auto systemMat = sparseXtDisc->getSystemMatrix();

    // sparse rhs assembly
    auto spatialNestedFEHierarchyForTemperature
            = sparseXtDisc->getSpatialNestedFEHierarchyForTemperature();
    auto spatialNestedFEHierarchyForHeatFlux
            = sparseXtDisc->getSpatialNestedFEHierarchyForHeatFlux();

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (temporalLevel, temporalLevel,
             spatialNestedFEHierarchyForTemperature,
             spatialNestedFEHierarchyForHeatFlux);

    auto dataSize = solutionHandler->getDataSize();
    auto dataOffsets = evalBlockOffsets(dataSize);
    auto rhs = std::make_shared<BlockVector>(dataOffsets);
    *rhs = 0.;
    sparseXtDisc->assembleRhs(rhs);

    // full system assembly
    std::unique_ptr<heat::LsqXtFEM> fullXtDisc
            = std::make_unique<heat::LsqXtFemH1Hdiv>
            (heatConfig, testCase);
    fullXtDisc->set(temporalMesh, spatialMesh);
    fullXtDisc->assembleSystem();
    fullXtDisc->buildSystemMatrix();
    auto trueSystemBlock11 = fullXtDisc->getSystemBlock11();
    auto trueSystemBlock12 = fullXtDisc->getSystemBlock12();
    auto trueSystemBlock21 = fullXtDisc->getSystemBlock21();
    auto trueSystemBlock22 = fullXtDisc->getSystemBlock22();
    auto trueSystemMat = fullXtDisc->getHeatMat();

    // full rhs assembly
    auto trueRhs = std::make_shared<BlockVector>(dataOffsets);
    *trueRhs = 0.;
    fullXtDisc->assembleRhs(trueRhs.get());

    // check
    double tol = 1E-8;

    trueSystemBlock11->Add(-1, *systemBlock11);
    ASSERT_LE(trueSystemBlock11->MaxNorm(), tol);

    trueSystemBlock12->Add(-1, *systemBlock12);
    ASSERT_LE(trueSystemBlock12->MaxNorm(), tol);

    trueSystemBlock21->Add(-1, *systemBlock21);
    ASSERT_LE(trueSystemBlock21->MaxNorm(), tol);

    trueSystemBlock22->Add(-1, *systemBlock22);
    ASSERT_LE(trueSystemBlock22->MaxNorm(), tol);

    trueSystemMat->Add(-1, *systemMat);
    ASSERT_LE(trueSystemMat->MaxNorm(), tol);

    trueRhs->Add(-1, *rhs);
    ASSERT_LE(trueRhs->Normlinf(), tol);
}


/**
 * @brief Tests the least-squares sparse space-time H1-Hdiv discretisation with unitSquare_test3
 *
 * Compares the sparse and full discretisations for a single level, with test case unitSquare_test3
 */
TEST(SparseDiscretisation, singleLevelH1HdivUnitSquareTest3)
{
    // config files for full and sparse versions
    std::string heatConfigFile
            = "../config_files/unit_tests/"
              "sparse_heat_discretisation/heat_square_test3.json";
    std::string sparseHeatConfigFile
            = "../config_files/unit_tests/"
              "sparse_heat_discretisation/sparseHeat_square_test3.json";
    auto heatConfig = getGlobalConfig(heatConfigFile);
    auto sparseHeatConfig = getGlobalConfig(sparseHeatConfigFile);

    auto testCase = heat::makeTestCase(heatConfig);

    int deg = heatConfig["deg"];
    ASSERT_EQ(deg, 1);

    int numLevels = 1;
    int spatialLevel = heatConfig["spatial_level"];
    int temporalLevel = heatConfig["temporal_level"];
    ASSERT_EQ(numLevels, sparseHeatConfig["num_levels"]);
    ASSERT_EQ(spatialLevel, sparseHeatConfig["min_spatial_level"]);
    ASSERT_EQ(temporalLevel, sparseHeatConfig["min_temporal_level"]);

    // spatial mesh hierarchy
    std::string inputDir
            = "../tests/input/sparse_heat_discretisation/";
    const std::string spatialMeshFile = inputDir+"mesh_l0.mesh";
    auto spatialMesh
            = std::make_shared<Mesh>(spatialMeshFile.c_str());
    for (int m=0; m<spatialLevel; m++) {
        spatialMesh->UniformRefinement();
    }
    auto spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();
    spatialMeshHierarchy->addMesh(spatialMesh);
    spatialMeshHierarchy->finalize();

    // temporal mesh
    double endTime = heatConfig["end_time"];
    int Nt = static_cast<int>(std::pow(2, temporalLevel));
    auto temporalMesh = std::make_shared<Mesh>(Nt, endTime);

    // sparse system assembly
    std::unique_ptr<sparseHeat::LsqSparseXtFem> sparseXtDisc
            = std::make_unique<sparseHeat::LsqSparseXtFemH1Hdiv>
            (sparseHeatConfig, testCase,
             numLevels, temporalLevel,
             spatialMeshHierarchy);
    sparseXtDisc->assembleSystemSubMatrices();
    sparseXtDisc->buildSystemMatrix();
    auto systemBlock11 = sparseXtDisc->getSystemBlock11();
    auto systemBlock12 = sparseXtDisc->getSystemBlock12();
    auto systemBlock21 = sparseXtDisc->getSystemBlock21();
    auto systemBlock22 = sparseXtDisc->getSystemBlock22();
    auto systemMat = sparseXtDisc->getSystemMatrix();

    // sparse rhs assembly
    auto spatialNestedFEHierarchyForTemperature
            = sparseXtDisc->getSpatialNestedFEHierarchyForTemperature();
    auto spatialNestedFEHierarchyForHeatFlux
            = sparseXtDisc->getSpatialNestedFEHierarchyForHeatFlux();

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (temporalLevel, temporalLevel,
             spatialNestedFEHierarchyForTemperature,
             spatialNestedFEHierarchyForHeatFlux);

    auto dataSize = solutionHandler->getDataSize();
    auto dataOffsets = evalBlockOffsets(dataSize);
    auto rhs = std::make_shared<BlockVector>(dataOffsets);
    *rhs = 0.;
    sparseXtDisc->assembleRhs(rhs);

    // full system assembly
    std::unique_ptr<heat::LsqXtFEM> fullXtDisc
            = std::make_unique<heat::LsqXtFemH1Hdiv>
            (heatConfig, testCase);
    fullXtDisc->set(temporalMesh, spatialMesh);
    fullXtDisc->assembleSystem();
    fullXtDisc->buildSystemMatrix();
    auto trueSystemBlock11 = fullXtDisc->getSystemBlock11();
    auto trueSystemBlock12 = fullXtDisc->getSystemBlock12();
    auto trueSystemBlock21 = fullXtDisc->getSystemBlock21();
    auto trueSystemBlock22 = fullXtDisc->getSystemBlock22();
    auto trueSystemMat = fullXtDisc->getHeatMat();

    // full rhs assembly
    auto trueRhs = std::make_shared<BlockVector>(dataOffsets);
    *trueRhs = 0.;
    fullXtDisc->assembleRhs(trueRhs.get());

    // check
    double tol = 1E-8;

    trueSystemBlock11->Add(-1, *systemBlock11);
    ASSERT_LE(trueSystemBlock11->MaxNorm(), tol);

    trueSystemBlock12->Add(-1, *systemBlock12);
    ASSERT_LE(trueSystemBlock12->MaxNorm(), tol);

    trueSystemBlock21->Add(-1, *systemBlock21);
    ASSERT_LE(trueSystemBlock21->MaxNorm(), tol);

    trueSystemBlock22->Add(-1, *systemBlock22);
    ASSERT_LE(trueSystemBlock22->MaxNorm(), tol);

    trueSystemMat->Add(-1, *systemMat);
    ASSERT_LE(trueSystemMat->MaxNorm(), tol);

    trueRhs->Add(-1, *rhs);
    ASSERT_LE(trueRhs->Normlinf(), tol);
}


/**
 * @brief Dry run of the least-squares sparse space-time H1-H1 discretisation
 */
TEST(SparseDiscretisation, buildH1H1)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/sparse_heat_discretisation/sparseHeat_dummy.json";
    auto config = getGlobalConfig(configFile);
    auto testCase = heat::makeTestCase(config);

    int deg = config["deg"];
    ASSERT_EQ(deg, 1);

    // spatial mesh hierarchy
    std::string input_dir
            = "../tests/input/sparse_heat_assembly/";
    const std::string spatialMeshFile1 = input_dir+"mesh_lx0";
    const std::string spatialMeshFile2 = input_dir+"mesh_lx1";
    const std::string spatialMeshFile3 = input_dir+"mesh_lx2";
    auto spatialMesh1 = std::make_shared<Mesh>(spatialMeshFile1.c_str());
    auto spatialMesh2 = std::make_shared<Mesh>(spatialMeshFile2.c_str());
    auto spatialMesh3 = std::make_shared<Mesh>(spatialMeshFile3.c_str());
    ASSERT_EQ(spatialMesh1->GetNE(), 1);
    ASSERT_EQ(spatialMesh2->GetNE(), 2);
    ASSERT_EQ(spatialMesh3->GetNE(), 3);

    auto spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();
    spatialMeshHierarchy->addMesh(spatialMesh1);
    spatialMeshHierarchy->addMesh(spatialMesh2);
    spatialMeshHierarchy->addMesh(spatialMesh3);
    spatialMeshHierarchy->finalize();

    // temporal hierarchy info
    int minTemporalLevel = 1;
    if (config.contains("min_temporal_level")) {
        minTemporalLevel = config["min_temporal_level"];
    }
    int numLevels = 1;
    if (config.contains("num_levels")) {
        numLevels = config["num_levels"];
    }
    int maxTemporalLevel = minTemporalLevel + numLevels - 1;
    ASSERT_EQ(minTemporalLevel, 1);
    ASSERT_EQ(numLevels, 3);
    ASSERT_EQ(maxTemporalLevel, 3);

    // discretisation
    std::unique_ptr<sparseHeat::LsqSparseXtFem> xtDisc
            = std::make_unique<sparseHeat::LsqSparseXtFemH1H1>
            (config, testCase, numLevels, minTemporalLevel);
    xtDisc->setNestedFEHierarchyAndSpatialBoundaryDofs
            (spatialMeshHierarchy);

    // system assembly
    xtDisc->assembleSystemSubMatrices();
    xtDisc->buildSystemMatrix();

    // rhs assembly
    auto spatialNestedFEHierarchyForTemperature
            = xtDisc->getSpatialNestedFEHierarchyForTemperature();
    auto spatialNestedFEHierarchyForHeatFlux
            = xtDisc->getSpatialNestedFEHierarchyForHeatFlux();

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (minTemporalLevel, maxTemporalLevel,
             spatialNestedFEHierarchyForTemperature,
             spatialNestedFEHierarchyForHeatFlux);

    auto dataSize = solutionHandler->getDataSize();
    auto dataOffsets = evalBlockOffsets(dataSize);
    auto rhs = std::make_shared<BlockVector>(dataOffsets);
    (*rhs) = 0;
    xtDisc->assembleRhs(rhs);

    // verify sizes
    auto systemBlock11 = xtDisc->getSystemBlock11();
    ASSERT_EQ(systemBlock11->NumRows(), 35);
    ASSERT_EQ(systemBlock11->NumCols(), 35);

    auto systemBlock12 = xtDisc->getSystemBlock12();
    ASSERT_EQ(systemBlock12->NumRows(), 35);
    ASSERT_EQ(systemBlock12->NumCols(), 70);

    auto systemBlock21 = xtDisc->getSystemBlock21();
    ASSERT_EQ(systemBlock21->NumRows(), 70);
    ASSERT_EQ(systemBlock21->NumCols(), 35);

    auto systemBlock22 = xtDisc->getSystemBlock22();
    ASSERT_EQ(systemBlock22->NumRows(), 70);
    ASSERT_EQ(systemBlock22->NumCols(), 70);

    auto systemMatrix = xtDisc->getSystemMatrix();
    ASSERT_EQ(systemMatrix->NumRows(), 105);
    ASSERT_EQ(systemMatrix->NumCols(), 105);
}

/**
 * @brief Tests the least-squares sparse space-time H1-H1 discretisation with unitSquare_test1
 *
 * Compares the sparse and full discretisations for a single level, with test case unitSquare_test1
 */
TEST(SparseDiscretisation, singleLevelH1H1UnitSquareTest1)
{
    // config files for full and sparse versions
    std::string heatConfigFile
            = "../config_files/unit_tests/"
              "sparse_heat_discretisation/heat_square_test1.json";
    std::string sparseHeatConfigFile
            = "../config_files/unit_tests/"
              "sparse_heat_discretisation/sparseHeat_square_test1.json";
    auto heatConfig = getGlobalConfig(heatConfigFile);
    auto sparseHeatConfig = getGlobalConfig(sparseHeatConfigFile);

    auto testCase = heat::makeTestCase(heatConfig);

    int deg = heatConfig["deg"];
    ASSERT_EQ(deg, 1);

    int numLevels = 1;
    int spatialLevel = heatConfig["spatial_level"];
    int temporalLevel = heatConfig["temporal_level"];
    ASSERT_EQ(numLevels, sparseHeatConfig["num_levels"]);
    ASSERT_EQ(spatialLevel, sparseHeatConfig["min_spatial_level"]);
    ASSERT_EQ(temporalLevel, sparseHeatConfig["min_temporal_level"]);

    // spatial mesh hierarchy
    std::string inputDir
            = "../tests/input/sparse_heat_discretisation/";
    const std::string spatialMeshFile = inputDir+"mesh_l0.mesh";
    auto spatialMesh
            = std::make_shared<Mesh>(spatialMeshFile.c_str());
    for (int m=0; m<spatialLevel; m++) {
        spatialMesh->UniformRefinement();
    }
    auto spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();
    spatialMeshHierarchy->addMesh(spatialMesh);
    spatialMeshHierarchy->finalize();

    // temporal mesh
    double endTime = heatConfig["end_time"];
    int Nt = static_cast<int>(std::pow(2, temporalLevel));
    auto temporalMesh = std::make_shared<Mesh>(Nt, endTime);

    // sparse system assembly
    std::unique_ptr<sparseHeat::LsqSparseXtFem> sparseXtDisc
            = std::make_unique<sparseHeat::LsqSparseXtFemH1H1>
            (sparseHeatConfig, testCase,
             numLevels, temporalLevel,
             spatialMeshHierarchy);
    sparseXtDisc->assembleSystemSubMatrices();
    sparseXtDisc->buildSystemMatrix();
    auto systemBlock11 = sparseXtDisc->getSystemBlock11();
    auto systemBlock12 = sparseXtDisc->getSystemBlock12();
    auto systemBlock21 = sparseXtDisc->getSystemBlock21();
    auto systemBlock22 = sparseXtDisc->getSystemBlock22();
    auto systemMat = sparseXtDisc->getSystemMatrix();

    // sparse rhs assembly
    auto spatialNestedFEHierarchyForTemperature
            = sparseXtDisc->getSpatialNestedFEHierarchyForTemperature();
    auto spatialNestedFEHierarchyForHeatFlux
            = sparseXtDisc->getSpatialNestedFEHierarchyForHeatFlux();

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (temporalLevel, temporalLevel,
             spatialNestedFEHierarchyForTemperature,
             spatialNestedFEHierarchyForHeatFlux);

    auto dataSize = solutionHandler->getDataSize();
    auto dataOffsets = evalBlockOffsets(dataSize);
    auto rhs = std::make_shared<BlockVector>(dataOffsets);
    *rhs = 0.;
    sparseXtDisc->assembleRhs(rhs);

    // full system assembly
    std::unique_ptr<heat::LsqXtFEM> fullXtDisc
            = std::make_unique<heat::LsqXtFemH1H1>
            (heatConfig, testCase);
    fullXtDisc->set(temporalMesh, spatialMesh);
    fullXtDisc->assembleSystem();
    fullXtDisc->buildSystemMatrix();
    auto trueSystemBlock11 = fullXtDisc->getSystemBlock11();
    auto trueSystemBlock12 = fullXtDisc->getSystemBlock12();
    auto trueSystemBlock21 = fullXtDisc->getSystemBlock21();
    auto trueSystemBlock22 = fullXtDisc->getSystemBlock22();
    auto trueSystemMat = fullXtDisc->getHeatMat();

    // full rhs assembly
    auto trueRhs = std::make_shared<BlockVector>(dataOffsets);
    *trueRhs = 0.;
    fullXtDisc->assembleRhs(trueRhs.get());

    // check
    double tol = 1E-8;

    trueSystemBlock11->Add(-1, *systemBlock11);
    ASSERT_LE(trueSystemBlock11->MaxNorm(), tol);

    trueSystemBlock12->Add(-1, *systemBlock12);
    ASSERT_LE(trueSystemBlock12->MaxNorm(), tol);

    trueSystemBlock21->Add(-1, *systemBlock21);
    ASSERT_LE(trueSystemBlock21->MaxNorm(), tol);

    trueSystemBlock22->Add(-1, *systemBlock22);
    ASSERT_LE(trueSystemBlock22->MaxNorm(), tol);

    trueSystemMat->Add(-1, *systemMat);
    ASSERT_LE(trueSystemMat->MaxNorm(), tol);

    trueRhs->Add(-1, *rhs);
    ASSERT_LE(trueRhs->Normlinf(), tol);
}
