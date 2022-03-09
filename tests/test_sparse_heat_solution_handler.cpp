#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../src/core/config.hpp"
#include "../src/sparse_heat/solution_handler.hpp"
#include "../src/sparse_heat/utilities.hpp"


using namespace mymfem;


/**
 * @brief Tests the class SolutionHandler for 2d mesh hierarchy
 */
TEST(SolutionHandler, dataSizes2d)
{
    std::string input_dir
            = "../tests/input/nested_hierarchy/2d/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 2);
    ASSERT_EQ(mesh2->GetNE(), 6);
    ASSERT_EQ(mesh3->GetNE(), 9);

    std::shared_ptr<NestedMeshHierarchy> meshHierarchy
            = std::make_shared<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);
    meshHierarchy->finalize();

    int dim = mesh1->Dimension();
    
    // H1 spaces for temperature
    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollH1);
    auto fes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollH1);
    auto fes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollH1);
    
    // RT spaces for heat-flux                                                   
    FiniteElementCollection *feCollRT
            = new RT_FECollection(0, dim);
    auto fes1RT = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollRT);
    auto fes2RT = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollRT);
    auto fes3RT = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollRT);

    // hierarchy of spatial FE spaces for temperature
    std::shared_ptr<NestedFEHierarchy> spatialNestedFEHierarchyTemperature
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);
    spatialNestedFEHierarchyTemperature->addFESpace(fes1H1);
    spatialNestedFEHierarchyTemperature->addFESpace(fes2H1);
    spatialNestedFEHierarchyTemperature->addFESpace(fes3H1);

    // hierarchy of spatial FE spaces for heat-flux
    std::shared_ptr<NestedFEHierarchy> spatialNestedFEHierarchyHeatFlux
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);
    spatialNestedFEHierarchyHeatFlux->addFESpace(fes1RT);
    spatialNestedFEHierarchyHeatFlux->addFESpace(fes2RT);
    spatialNestedFEHierarchyHeatFlux->addFESpace(fes3RT);

    int minLevel = 2;
    int maxLevel = 4;
    int numLevels = maxLevel - minLevel + 1;
    assert(numLevels == 3);

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (minLevel, maxLevel,
             spatialNestedFEHierarchyTemperature,
             spatialNestedFEHierarchyHeatFlux);

    auto& temperatureData = solutionHandler->getTemperatureData();
    auto& heatFluxData = solutionHandler->getHeatFluxData();
    auto temporalHierarchicalFESizes
            = evalTemporalBlockSizes(minLevel, maxLevel);

    int trueTemperatureDataSize
            = temporalHierarchicalFESizes[0]*fes3H1->GetTrueVSize()
            + temporalHierarchicalFESizes[1]*fes2H1->GetTrueVSize()
            + temporalHierarchicalFESizes[2]*fes1H1->GetTrueVSize();
    ASSERT_EQ(temperatureData.Size(), trueTemperatureDataSize);

    int trueHeatFluxDataSize
            = temporalHierarchicalFESizes[0]*fes3RT->GetTrueVSize()
            + temporalHierarchicalFESizes[1]*fes2RT->GetTrueVSize()
            + temporalHierarchicalFESizes[2]*fes1RT->GetTrueVSize();
    ASSERT_EQ(heatFluxData.Size(), trueHeatFluxDataSize);

    auto solutionData = solutionHandler->getData();
    ASSERT_EQ((*solutionData).Size(),
              trueTemperatureDataSize + trueHeatFluxDataSize);
}

/**
 * @brief Tests the class SolutionHandler for 23 mesh hierarchy
 */
TEST(SolutionHandler, dataSizes3d)
{
    std::string input_dir
            = "../tests/input/nested_hierarchy/3d/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 1);
    ASSERT_EQ(mesh2->GetNE(), 2);
    ASSERT_EQ(mesh3->GetNE(), 3);

    std::shared_ptr<NestedMeshHierarchy> meshHierarchy
            = std::make_shared<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);
    meshHierarchy->finalize();

    int dim = mesh1->Dimension();

    // H1 spaces for temperature
    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollH1);
    auto fes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollH1);
    auto fes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollH1);

    // RT spaces for heat-flux
    FiniteElementCollection *feCollRT
            = new RT_FECollection(0, dim);
    auto fes1RT = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollRT);
    auto fes2RT = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollRT);
    auto fes3RT = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollRT);

    // hierarchy of spatial FE spaces for temperature
    std::shared_ptr<NestedFEHierarchy> spatialNestedFEHierarchyTemperature
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);
    spatialNestedFEHierarchyTemperature->addFESpace(fes1H1);
    spatialNestedFEHierarchyTemperature->addFESpace(fes2H1);
    spatialNestedFEHierarchyTemperature->addFESpace(fes3H1);

    // hierarchy of spatial FE spaces for heat-flux
    std::shared_ptr<NestedFEHierarchy> spatialNestedFEHierarchyHeatFlux
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);
    spatialNestedFEHierarchyHeatFlux->addFESpace(fes1RT);
    spatialNestedFEHierarchyHeatFlux->addFESpace(fes2RT);
    spatialNestedFEHierarchyHeatFlux->addFESpace(fes3RT);

    int minLevel = 2;
    int maxLevel = 4;
    int numLevels = maxLevel - minLevel + 1;
    assert(numLevels == 3);

    auto solutionHandler
            = std::make_unique<sparseHeat::SolutionHandler>
            (minLevel, maxLevel,
             spatialNestedFEHierarchyTemperature,
             spatialNestedFEHierarchyHeatFlux);

    auto& temperatureData = solutionHandler->getTemperatureData();
    auto& heatFluxData = solutionHandler->getHeatFluxData();
    auto temporalHierarchicalFESizes
            = evalTemporalBlockSizes(minLevel, maxLevel);

    int trueTemperatureDataSize
            = temporalHierarchicalFESizes[0]*fes3H1->GetTrueVSize()
            + temporalHierarchicalFESizes[1]*fes2H1->GetTrueVSize()
            + temporalHierarchicalFESizes[2]*fes1H1->GetTrueVSize();
    ASSERT_EQ(temperatureData.Size(), trueTemperatureDataSize);

    int trueHeatFluxDataSize
            = temporalHierarchicalFESizes[0]*fes3RT->GetTrueVSize()
            + temporalHierarchicalFESizes[1]*fes2RT->GetTrueVSize()
            + temporalHierarchicalFESizes[2]*fes1RT->GetTrueVSize();
    ASSERT_EQ(heatFluxData.Size(), trueHeatFluxDataSize);

    auto solutionData = solutionHandler->getData();
    ASSERT_EQ((*solutionData).Size(),
              trueTemperatureDataSize + trueHeatFluxDataSize);
}
