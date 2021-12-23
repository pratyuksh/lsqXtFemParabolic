#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../src/core/config.hpp"
#include "../src/heat/assembly.hpp"
#include "../src/sparse_heat/spatial_assembly.hpp"


using namespace mymfem;


/**
 * @brief Tests the assembly of sparse Vector mass matrix
 */
TEST(SparseSpatialAssembly, vectorMassMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_heat_assembly/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 1);
    ASSERT_EQ(mesh2->GetNE(), 2);
    ASSERT_EQ(mesh3->GetNE(), 3);

    auto meshHierarchy
            = std::make_shared<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);
    meshHierarchy->finalize();

    auto nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl, dim);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl, dim);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl, dim);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto numDims = nestedFEHierarchy->getNumDims();
    ASSERT_EQ(numDims[0], 6);
    ASSERT_EQ(numDims[1], 8);
    ASSERT_EQ(numDims[2], 10);

    auto massBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> massIntegator
            = std::make_shared<sparseHeat::SpatialVectorMassIntegrator>();

    massBilinearForm->addDomainIntegrator(massIntegator);
    massBilinearForm->assemble();

    auto blockMassMatrix
            = massBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        BilinearForm massBf(fes1.get());
        massBf.AddDomainIntegrator(new VectorMassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass00 = massBf.LoseMat();
        trueBlockMass00->Add(-1, blockMassMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockMass00->MaxNorm(), tol);
    }
    {
        BilinearForm massBf(fes2.get());
        massBf.AddDomainIntegrator(new VectorMassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass11 = massBf.LoseMat();
        trueBlockMass11->Add(-1, blockMassMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockMass11->MaxNorm(), tol);
    }
    {
        BilinearForm massBf(fes3.get());
        massBf.AddDomainIntegrator(new VectorMassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass22 = massBf.LoseMat();
        trueBlockMass22->Add(-1, blockMassMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockMass22->MaxNorm(), tol);
    }

//    for(int i=1; i<blockMassMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<i; j++) {
//            blockMassMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    delete feColl;
}


/**
 * @brief Tests the assembly of sparse Vector stiffness matrix
 */
TEST(SparseSpatialAssembly, vectorStiffnessMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_heat_assembly/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 1);
    ASSERT_EQ(mesh2->GetNE(), 2);
    ASSERT_EQ(mesh3->GetNE(), 3);

    auto meshHierarchy
            = std::make_shared<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);
    meshHierarchy->finalize();

    auto nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl, dim);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl, dim);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl, dim);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto numDims = nestedFEHierarchy->getNumDims();
    ASSERT_EQ(numDims[0], 6);
    ASSERT_EQ(numDims[1], 8);
    ASSERT_EQ(numDims[2], 10);

    auto stiffnessBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> stiffnessIntegator
            = std::make_shared<sparseHeat::SpatialVectorStiffnessIntegrator>();

    stiffnessBilinearForm->addDomainIntegrator(stiffnessIntegator);
    stiffnessBilinearForm->assemble();

    auto blockStiffnessMatrix = stiffnessBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        BilinearForm stiffnessBf(fes1.get());
        stiffnessBf.AddDomainIntegrator
                (new heat::SpatialVectorStiffnessIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness00 = stiffnessBf.LoseMat();
        trueBlockStiffness00->
                Add(-1, blockStiffnessMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockStiffness00->MaxNorm(), tol);
    }
    {
        BilinearForm stiffnessBf(fes2.get());
        stiffnessBf.AddDomainIntegrator
                (new heat::SpatialVectorStiffnessIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness11 = stiffnessBf.LoseMat();
        trueBlockStiffness11->
                Add(-1, blockStiffnessMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockStiffness11->MaxNorm(), tol);
    }
    {
        BilinearForm stiffnessBf(fes3.get());
        stiffnessBf.AddDomainIntegrator
                (new heat::SpatialVectorStiffnessIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness22 = stiffnessBf.LoseMat();
        trueBlockStiffness22->
                Add(-1, blockStiffnessMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockStiffness22->MaxNorm(), tol);
    }

//    for(int i=1; i<blockStiffnessMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<i; j++) {
//            blockStiffnessMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    delete feColl;
}


/**
 * @brief Tests the assembly of sparse Vector gradient matrix
 */
TEST(SparseSpatialAssembly, vectorGradientMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_heat_assembly/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 1);
    ASSERT_EQ(mesh2->GetNE(), 2);
    ASSERT_EQ(mesh3->GetNE(), 3);

    auto meshHierarchy
            = std::make_shared<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);
    meshHierarchy->finalize();

    int dim = mesh1->Dimension();

    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);

    // Vector H1 spaces
    auto vfes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                        feCollH1, dim);
    auto vfes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                        feCollH1, dim);
    auto vfes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                        feCollH1, dim);

    // H1 spaces
    auto fes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollH1);
    auto fes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollH1);
    auto fes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollH1);

    // hierarchy of test spaces
    auto testNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);
    testNestedFEHierarchy->addFESpace(vfes1H1);
    testNestedFEHierarchy->addFESpace(vfes2H1);
    testNestedFEHierarchy->addFESpace(vfes3H1);

    auto testNumDims = testNestedFEHierarchy->getNumDims();
    ASSERT_EQ(testNumDims[0], 6);
    ASSERT_EQ(testNumDims[1], 8);
    ASSERT_EQ(testNumDims[2], 10);

    // hierarchy of trial spaces
    auto trialNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);
    trialNestedFEHierarchy->addFESpace(fes1H1);
    trialNestedFEHierarchy->addFESpace(fes2H1);
    trialNestedFEHierarchy->addFESpace(fes3H1);

    // assemble gradient mixed bilinear form
    auto gradientBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (trialNestedFEHierarchy, testNestedFEHierarchy);

    std::shared_ptr<BlockMixedBilinearFormIntegrator> gradientIntegator
            = std::make_shared<sparseHeat::SpatialVectorGradientIntegrator>();

    gradientBilinearForm->addDomainIntegrator(gradientIntegator);
    gradientBilinearForm->assemble();

    auto blockGradientMatrix
            = gradientBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        MixedBilinearForm gradientBf(fes1H1.get(),
                                     vfes1H1.get());
        gradientBf.AddDomainIntegrator
                (new heat::SpatialVectorGradientIntegrator{});
        gradientBf.Assemble();
        auto trueBlockGradient00 = gradientBf.LoseMat();
        trueBlockGradient00->
                Add(-1, blockGradientMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockGradient00->MaxNorm(), tol);
    }
    {
        MixedBilinearForm gradientBf(fes2H1.get(),
                                     vfes2H1.get());
        gradientBf.AddDomainIntegrator
                (new heat::SpatialVectorGradientIntegrator{});
        gradientBf.Assemble();
        auto trueBlockGradient11 = gradientBf.LoseMat();
        trueBlockGradient11->
                Add(-1, blockGradientMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockGradient11->MaxNorm(), tol);
    }
    {
        MixedBilinearForm gradientBf(fes3H1.get(),
                                     vfes3H1.get());
        gradientBf.AddDomainIntegrator
                (new heat::SpatialVectorGradientIntegrator{});
        gradientBf.Assemble();
        auto trueBlockGradient22 = gradientBf.LoseMat();
        trueBlockGradient22->
                Add(-1, blockGradientMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockGradient22->MaxNorm(), tol);
    }

//    for(int i=1; i<blockGradientMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<i; j++) {
//            blockGradientMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    delete feCollH1;
}


/**
 * @brief Tests the assembly of sparse Vector divergence matrix
 */
TEST(SparseSpatialAssembly, vectorDivergenceMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_heat_assembly/";

    const std::string meshFile1 = input_dir+"mesh_lx0";
    const std::string meshFile2 = input_dir+"mesh_lx1";
    const std::string meshFile3 = input_dir+"mesh_lx2";

    auto mesh1 = std::make_shared<Mesh>(meshFile1.c_str());
    auto mesh2 = std::make_shared<Mesh>(meshFile2.c_str());
    auto mesh3 = std::make_shared<Mesh>(meshFile3.c_str());

    ASSERT_EQ(mesh1->GetNE(), 1);
    ASSERT_EQ(mesh2->GetNE(), 2);
    ASSERT_EQ(mesh3->GetNE(), 3);

    auto meshHierarchy
            = std::make_shared<NestedMeshHierarchy>();
    meshHierarchy->addMesh(mesh1);
    meshHierarchy->addMesh(mesh2);
    meshHierarchy->addMesh(mesh3);
    meshHierarchy->finalize();

    int dim = mesh1->Dimension();

    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);

    // Vector H1 spaces
    auto vfes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                        feCollH1,
                                                        dim);
    auto vfes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                        feCollH1,
                                                        dim);
    auto vfes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                        feCollH1,
                                                        dim);

    // H1 spaces
    auto fes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollH1);
    auto fes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollH1);
    auto fes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollH1);

    // hierarchy of test spaces
    auto testNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);
    testNestedFEHierarchy->addFESpace(fes1H1);
    testNestedFEHierarchy->addFESpace(fes2H1);
    testNestedFEHierarchy->addFESpace(fes3H1);

    // hierarchy of trial spaces
    auto trialNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);
    trialNestedFEHierarchy->addFESpace(vfes1H1);
    trialNestedFEHierarchy->addFESpace(vfes2H1);
    trialNestedFEHierarchy->addFESpace(vfes3H1);

    auto trialNumDims = trialNestedFEHierarchy->getNumDims();
    ASSERT_EQ(trialNumDims[0], 6);
    ASSERT_EQ(trialNumDims[1], 8);
    ASSERT_EQ(trialNumDims[2], 10);

    // assemble divergence mixed bilinear form
    auto divergenceBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (trialNestedFEHierarchy, testNestedFEHierarchy);

    std::shared_ptr<BlockMixedBilinearFormIntegrator> divergenceIntegator
            = std::make_shared<sparseHeat::SpatialVectorDivergenceIntegrator>();

    divergenceBilinearForm->addDomainIntegrator(divergenceIntegator);
    divergenceBilinearForm->assemble();

    auto blockDivergenceMatrix
            = divergenceBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        MixedBilinearForm divergenceBf(vfes1H1.get(),
                                       fes1H1.get());
        divergenceBf.AddDomainIntegrator
                (new VectorDivergenceIntegrator{});
        divergenceBf.Assemble();
        auto trueBlockDivergence00 = divergenceBf.LoseMat();
        trueBlockDivergence00->
                Add(-1, blockDivergenceMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockDivergence00->MaxNorm(), tol);
    }
    {
        MixedBilinearForm divergenceBf(vfes2H1.get(),
                                       fes2H1.get());
        divergenceBf.AddDomainIntegrator
                (new VectorDivergenceIntegrator{});
        divergenceBf.Assemble();
        auto trueBlockDivergence11 = divergenceBf.LoseMat();
        trueBlockDivergence11->
                Add(-1, blockDivergenceMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockDivergence11->MaxNorm(), tol);
    }
    {
        MixedBilinearForm divergenceBf(vfes3H1.get(),
                                       fes3H1.get());
        divergenceBf.AddDomainIntegrator
                (new VectorDivergenceIntegrator{});
        divergenceBf.Assemble();
        auto trueBlockDivergence22 = divergenceBf.LoseMat();
        trueBlockDivergence22->
                Add(-1, blockDivergenceMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockDivergence22->MaxNorm(), tol);
    }

//    for(int i=1; i<blockDivergenceMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<i; j++) {
//            blockDivergenceMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    delete feCollH1;
}
