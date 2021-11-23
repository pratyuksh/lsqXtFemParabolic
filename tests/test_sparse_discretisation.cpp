#include <gtest/gtest.h>

#include "mfem.hpp"
#include <iostream>

#include "../src/heat/test_cases_factory.hpp"
#include "../src/sparse_heat/discretisation.hpp"
#include "../src/sparse_heat/solution_handler.hpp"
#include "../src/sparse_heat/utilities.hpp"

using namespace mfem;

/**
 * @brief Tests the least-squares space-time H1-Hdiv discretisation
 */
TEST(SparseDiscretisation, buildH1Hdiv)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/sparseHeat_dummy.json";
    auto config = getGlobalConfig(configFile);
    auto testCase = heat::makeTestCase(config);

    int deg = config["deg"];
    ASSERT_EQ(deg, 1);

    // spatial mesh hierarchy
    std::string input_dir
            = "../tests/input/sparse_assembly/";
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
            (config, testCase);
    xtDisc->setNestedFEHierarchyAndSpatialBoundaryDofs
            (spatialMeshHierarchy);

    // system assembly
    xtDisc->assembleSystemSubMatrices();
    xtDisc->buildSystemMatrix();

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

    //systemBlock11->Print();

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
    auto B = std::make_shared<BlockVector>(dataOffsets);
    xtDisc->assembleRhs(B);
}
