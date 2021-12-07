#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>
#include <chrono>

#include "../src/core/config.hpp"
#include "../src/mymfem/nested_hierarchy.hpp"
#include "../src/mymfem/my_bilinearForms.hpp"


using namespace mymfem;

/**
 * @brief Tests the class BlockBilinearForm
 */
TEST(MyBilinearForms, BlockBilinearForm)
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

    std::shared_ptr<NestedFEHierarchy> nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    std::unique_ptr<BlockBilinearForm> blockBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    Array<int> trueBlockOffsets(4);
    trueBlockOffsets[0] = 0;
    trueBlockOffsets[1] = fes1->GetTrueVSize();
    trueBlockOffsets[2] = fes2->GetTrueVSize() + trueBlockOffsets[1];
    trueBlockOffsets[3] = fes3->GetTrueVSize() + trueBlockOffsets[2];

    // check block offsets
    auto blockOffsets = blockBilinearForm->getBlockOffsets();
    ASSERT_EQ(blockOffsets, trueBlockOffsets);

    // check that memory has been allocated for the block matrix
    blockBilinearForm->allocateBlockMatrix();
    auto blockMatrix = blockBilinearForm->getBlockMatrix();

    ASSERT_EQ(blockMatrix->NumRowBlocks(), 3);
    ASSERT_EQ(blockMatrix->NumColBlocks(), 3);
    for(int i=0; i<blockMatrix->NumRowBlocks(); i++)
        for (int j=0; j<blockMatrix->NumColBlocks(); j++) {
            ASSERT_NE(&blockMatrix->GetBlock(i,j), nullptr);
        }
}


/**
 * @brief Tests the class BlockMixedBilinearForm
 */
TEST(MyBilinearForms, BlockMixedBilinearForm)
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

    // hierarchy of test spaces
    std::shared_ptr<NestedFEHierarchy> testNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);
    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollH1);
    auto fes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollH1);
    auto fes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollH1);
    testNestedFEHierarchy->addFESpace(fes1H1);
    testNestedFEHierarchy->addFESpace(fes2H1);
    testNestedFEHierarchy->addFESpace(fes3H1);

    // hierarchy of trial spaces
    std::shared_ptr<NestedFEHierarchy> trialNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy> (meshHierarchy);
    FiniteElementCollection *feCollRT
            = new RT_FECollection(0, dim);
    auto fes1RT = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollRT);
    auto fes2RT = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollRT);
    auto fes3RT = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollRT);
    trialNestedFEHierarchy->addFESpace(fes1RT);
    trialNestedFEHierarchy->addFESpace(fes2RT);
    trialNestedFEHierarchy->addFESpace(fes3RT);

    std::unique_ptr<BlockMixedBilinearForm> blockMixedBilinearForm
            = std::make_unique<BlockMixedBilinearForm>
            (trialNestedFEHierarchy, testNestedFEHierarchy);

    // check row block offsets
    Array<int> trueRowBlockOffsets(4);
    trueRowBlockOffsets[0] = 0;
    trueRowBlockOffsets[1] = fes1H1->GetTrueVSize();
    trueRowBlockOffsets[2] = fes2H1->GetTrueVSize() + trueRowBlockOffsets[1];
    trueRowBlockOffsets[3] = fes3H1->GetTrueVSize() + trueRowBlockOffsets[2];

    auto rowBlockOffsets = blockMixedBilinearForm->getRowBlockOffsets();
    ASSERT_EQ(rowBlockOffsets, trueRowBlockOffsets);

    // check column block offsets
    Array<int> trueColBlockOffsets(4);
    trueColBlockOffsets[0] = 0;
    trueColBlockOffsets[1] = fes1RT->GetTrueVSize();
    trueColBlockOffsets[2] = fes2RT->GetTrueVSize() + trueColBlockOffsets[1];
    trueColBlockOffsets[3] = fes3RT->GetTrueVSize() + trueColBlockOffsets[2];

    auto colBlockOffsets = blockMixedBilinearForm->getColBlockOffsets();
    ASSERT_EQ(colBlockOffsets, trueColBlockOffsets);

    // check that memory has been allocated for the block matrix
    blockMixedBilinearForm->allocateBlockMatrix();
    auto blockMatrix = blockMixedBilinearForm->getBlockMatrix();

    ASSERT_EQ(blockMatrix->NumRowBlocks(), 3);
    ASSERT_EQ(blockMatrix->NumColBlocks(), 3);
    for(int i=0; i<blockMatrix->NumRowBlocks(); i++)
        for (int j=0; j<blockMatrix->NumColBlocks(); j++) {
            ASSERT_NE(&blockMatrix->GetBlock(i,j), nullptr);
        }
}
