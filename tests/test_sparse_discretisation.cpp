#include <gtest/gtest.h>

#include "mfem.hpp"
#include <iostream>

#include "../src/heat/test_cases_factory.hpp"
#include "../src/sparse_heat/discretisation.hpp"

using namespace mfem;

/**
 * @brief Tests the least-squares space-time H1-Hdiv discretisation
 */
TEST(SparseDiscretisation, H1Hdiv)
{
    // config and test case
    std::string configFile
            = "../config_files/unit_tests/sparseHeat_square_test1.json";
    auto config = getGlobalConfig(configFile);
    auto testCase = heat::makeTestCase(config);

    // spatial mesh hierarchy
    std::string input_dir
            = "../tests/input/sparse_assembly/";
    const std::string spatialMeshFile1 = input_dir+"mesh_lx0";
    const std::string spatialMeshFile2 = input_dir+"mesh_lx1";
    const std::string spatialMeshFile3 = input_dir+"mesh_lx2";
    auto spatialMesh1 = std::make_shared<Mesh>(spatialMeshFile1.c_str());
    auto spatialMesh2 = std::make_shared<Mesh>(spatialMeshFile2.c_str());
    auto spatialMesh3 = std::make_shared<Mesh>(spatialMeshFile3.c_str());
    auto spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();
    spatialMeshHierarchy->addMesh(spatialMesh1);
    spatialMeshHierarchy->addMesh(spatialMesh2);
    spatialMeshHierarchy->addMesh(spatialMesh3);
    spatialMeshHierarchy->buildHierarchicalTranformations();

    // discretisation
    std::unique_ptr<sparseHeat::LsqSparseXtFem> xtDisc
            = std::make_unique<sparseHeat::LsqSparseXtFemH1Hdiv>
            (config, testCase);

    xtDisc->setNestedFEHierarchy(spatialMeshHierarchy);
    xtDisc->assembleSystemSubMatrices();
    xtDisc->buildSystemBlocks();

    auto block11 = xtDisc->getSystemBlock11();
    //block11->Print();
}
