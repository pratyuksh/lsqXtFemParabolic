#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>
#include <chrono>

#include "../src/core/config.hpp"
#include "../src/mymfem/mesh_hierarchy.hpp"

using namespace mymfem;

/**
 * @brief Tests the construction of parent and child hierarchy
 * between two nested meshes.
 */
TEST(MeshHierarchy, singleLevel)
{
    std::string input_dir
            = "../input/mesh_hierarchy/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 2);
    ASSERT_EQ(mesh2->GetNE(), 6);
    ASSERT_EQ(mesh3->GetNE(), 9);

    // generate successive hierarchies
    MeshHierarchy meshHierarchy12(mesh1, mesh2);
    meshHierarchy12.build();
    auto hierarchy12 = meshHierarchy12.get();

    MeshHierarchy meshHierarchy23(mesh2, mesh3);
    meshHierarchy23.build();
    auto hierarchy23 = meshHierarchy23.get();

    auto trueHierarchy12
            = std::make_unique<MeshHierarchyTable>(mesh1->GetNE());
    // children of the element 0 in coarse mesh
    trueHierarchy12->add(0, 0);
    trueHierarchy12->add(0, 3);

    // children of the element 1 in coarse mesh
    trueHierarchy12->add(1, 4);
    trueHierarchy12->add(1, 2);
    trueHierarchy12->add(1, 1);
    trueHierarchy12->add(1, 5);

    for (int i=0; i<mesh1->GetNE(); i++) {
        auto buf = (*hierarchy12)(i);
        auto trueBuf = (*trueHierarchy12)(i);
        ASSERT_EQ(buf, trueBuf);
    }

    auto trueHierarchy23
            = std::make_unique<MeshHierarchyTable>(mesh2->GetNE());
    // children of the element 0 in coarse mesh
    trueHierarchy23->add(0, 0);
    trueHierarchy23->add(0, 1);
    trueHierarchy23->add(0, 2);

    // children of the element 1 in coarse mesh
    trueHierarchy23->add(1, 5);

    // children of the element 2 in coarse mesh
    trueHierarchy23->add(2, 6);

    // children of the element 3 in coarse mesh
    trueHierarchy23->add(3, 3);
    trueHierarchy23->add(3, 4);

    // children of the element 4 in coarse mesh
    trueHierarchy23->add(4, 7);

    // children of the element 5 in coarse mesh
    trueHierarchy23->add(5, 8);

    for (int i=0; i<mesh2->GetNE(); i++) {
        auto buf = (*hierarchy23)(i);
        auto trueBuf = (*trueHierarchy23)(i);
        ASSERT_EQ(buf, trueBuf);
    }
}

/**
 * @brief Tests the construction of parent and child hierarchy
 * across multiple levels
 */
TEST(MeshHierarchy, multiLevel)
{
    std::string input_dir
            = "../input/mesh_hierarchy/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 2);
    ASSERT_EQ(mesh2->GetNE(), 6);
    ASSERT_EQ(mesh3->GetNE(), 9);

    // generate successive hierarchies
    MeshHierarchy meshHierarchy12(mesh1, mesh2);
    meshHierarchy12.build();
    auto hierarchy12 = meshHierarchy12.get();

    MeshHierarchy meshHierarchy23(mesh2, mesh3);
    meshHierarchy23.build();
    auto hierarchy23 = meshHierarchy23.get();

    // collect hierarchies between successive levels
    std::vector<std::shared_ptr<MeshHierarchyTable>>
            hierarchies = {hierarchy12, hierarchy23};

    // generate hierarchy across multiple levels
    auto hierarchy13 = buildMLMeshHierarchy(hierarchies);

    // ground truth
    MeshHierarchy meshHierarchy13(mesh1, mesh3);
    meshHierarchy13.build();
    auto trueHierarchy13 = meshHierarchy13.get();

    for (int i=0; i<mesh1->GetNE(); i++) {
        auto buf = (*hierarchy13)(i);
        auto trueBuf = (*trueHierarchy13)(i);
        ASSERT_EQ(buf, trueBuf);
    }
    hierarchy13->print();
}
