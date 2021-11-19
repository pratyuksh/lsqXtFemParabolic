#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>
#include <chrono>

#include "../src/core/config.hpp"
#include "../src/mymfem/nested_hierarchy.hpp"

using namespace mymfem;

/**
 * @brief Tests the construction of parent and child hierarchy
 * between two nested meshes.
 */
TEST(NestedMeshHierarchy, singleLevel)
{
    std::string input_dir
            = "../tests/input/nested_hierarchy/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 2);
    ASSERT_EQ(mesh2->GetNE(), 6);
    ASSERT_EQ(mesh3->GetNE(), 9);

    std::unique_ptr<NestedMeshHierarchy> meshHierarchy
            = std::make_unique<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);

    meshHierarchy->buildHierarchicalTranformations();

    auto hierarchicalTransformations = meshHierarchy->getTransformations();

    auto trueHierarchicalTransformationBetweenMeshes12
            = std::make_unique<HierarchicalMeshTransformationTable>(mesh1->GetNE());
    // children of the element 0 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes12->add(0, 0);
    trueHierarchicalTransformationBetweenMeshes12->add(0, 3);

    // children of the element 1 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes12->add(1, 4);
    trueHierarchicalTransformationBetweenMeshes12->add(1, 2);
    trueHierarchicalTransformationBetweenMeshes12->add(1, 1);
    trueHierarchicalTransformationBetweenMeshes12->add(1, 5);

    for (int i=0; i<mesh1->GetNE(); i++) {
        auto buf = (*hierarchicalTransformations[0])(i);
        auto trueBuf = (*trueHierarchicalTransformationBetweenMeshes12)(i);
        ASSERT_EQ(buf, trueBuf);
    }

    auto trueHierarchicalTransformationBetweenMeshes23
            = std::make_unique<HierarchicalMeshTransformationTable>(mesh2->GetNE());
    // children of the element 0 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes23->add(0, 0);
    trueHierarchicalTransformationBetweenMeshes23->add(0, 1);
    trueHierarchicalTransformationBetweenMeshes23->add(0, 2);

    // children of the element 1 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes23->add(1, 5);

    // children of the element 2 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes23->add(2, 6);

    // children of the element 3 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes23->add(3, 3);
    trueHierarchicalTransformationBetweenMeshes23->add(3, 4);

    // children of the element 4 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes23->add(4, 7);

    // children of the element 5 in coarse mesh
    trueHierarchicalTransformationBetweenMeshes23->add(5, 8);

    for (int i=0; i<mesh2->GetNE(); i++) {
        auto buf = (*hierarchicalTransformations[1])(i);
        auto trueBuf = (*trueHierarchicalTransformationBetweenMeshes23)(i);
        ASSERT_EQ(buf, trueBuf);
    }
}

/**
 * @brief Tests the construction of parent and child hierarchy
 * across multiple levels
 */
TEST(NestedMeshHierarchy, multiLevel)
{
    std::string input_dir
            = "../tests/input/nested_hierarchy/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 2);
    ASSERT_EQ(mesh2->GetNE(), 6);
    ASSERT_EQ(mesh3->GetNE(), 9);

    std::unique_ptr<NestedMeshHierarchy> meshHierarchy
            = std::make_unique<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);

    meshHierarchy->buildHierarchicalTranformations();

    auto hierarchicalTransformations
            = meshHierarchy->getTransformations();

    // generate hierarchy across multiple levels

    // ground truth
    std::unique_ptr<NestedMeshHierarchy> trueMeshHierarchy
            = std::make_unique<NestedMeshHierarchy>();
    trueMeshHierarchy->addMesh(mesh1);
    trueMeshHierarchy->addMesh(mesh3);

    trueMeshHierarchy->buildHierarchicalTranformations();

    auto trueHierarchicalTransformations
            = trueMeshHierarchy->getTransformations();
    trueHierarchicalTransformations[0]->print();

    {
        auto hierarchicalTransformationMeshes13
                = buildMultiLevelMeshHierarchy(hierarchicalTransformations);

        for (int i=0; i<mesh1->GetNE(); i++) {
            auto buf = (*hierarchicalTransformationMeshes13)(i);
            auto trueBuf = (*trueHierarchicalTransformations[0])(i);
            ASSERT_EQ(buf, trueBuf);
        }
    }

    {
        auto hierarchicalTransformationMeshes13
                = buildMultiLevelMeshHierarchy
                (hierarchicalTransformations[0],
                hierarchicalTransformations[1]);

        for (int i=0; i<mesh1->GetNE(); i++) {
            auto buf = (*hierarchicalTransformationMeshes13)(i);
            auto trueBuf = (*trueHierarchicalTransformations[0])(i);
            ASSERT_EQ(buf, trueBuf);
        }
    }
}

/**
 * @brief Tests the construction of nested FE hierarchy
 */
TEST(NestedFEHierarchy, build)
{
    std::string input_dir
            = "../tests/input/nested_hierarchy/";

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

    meshHierarchy->buildHierarchicalTranformations();

    std::unique_ptr<NestedFEHierarchy> nestedFEHierarchy
            = std::make_unique<NestedFEHierarchy> (meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);

    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto numLevels = nestedFEHierarchy->getNumLevels();
    auto feSpaces = nestedFEHierarchy->getFESpaces();
    auto numFEDims = nestedFEHierarchy->getNumDims();

    ASSERT_EQ(numLevels, 3);

    ASSERT_EQ(feSpaces[0]->GetMesh(), mesh1.get());
    ASSERT_EQ(feSpaces[1]->GetMesh(), mesh2.get());
    ASSERT_EQ(feSpaces[2]->GetMesh(), mesh3.get());

    ASSERT_EQ(numFEDims[0], mesh1->GetNV());
    ASSERT_EQ(numFEDims[1], mesh2->GetNV());
    ASSERT_EQ(numFEDims[2], mesh3->GetNV());
}
