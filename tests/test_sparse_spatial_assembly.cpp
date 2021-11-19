#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../src/core/config.hpp"
#include "../src/heat/assembly.hpp"
#include "../src/sparse_heat/spatial_assembly.hpp"


using namespace mymfem;

/**
 * @brief Tests the assembly of sparse mass matrix
 */
TEST(SparseSpatialAssembly, massMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

    auto nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto massBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> massIntegator
            = std::make_shared<sparseHeat::SpatialMassIntegrator>();

    massBilinearForm->addDomainIntegrator(massIntegator);
    massBilinearForm->assemble();

    auto blockMassMatrix
            = massBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        BilinearForm massBf(fes1.get());
        massBf.AddDomainIntegrator(new MassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass00 = massBf.LoseMat();
        trueBlockMass00->Add(-1, blockMassMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockMass00->MaxNorm(), tol);
    }
    {
        BilinearForm massBf(fes2.get());
        massBf.AddDomainIntegrator(new MassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass11 = massBf.LoseMat();
        trueBlockMass11->Add(-1, blockMassMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockMass11->MaxNorm(), tol);
    }
    {
        BilinearForm massBf(fes3.get());
        massBf.AddDomainIntegrator(new MassIntegrator{});
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

    // check block(1,0)
    {
        auto mat = blockMassMatrix->GetBlock(1,0);

        ASSERT_LE(std::abs(mat.Elem(0,0) - 1./12), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - 1./24), tol);
        ASSERT_LE(std::abs(mat.Elem(0,2) - 1./24), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 5./96), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - 1./96), tol);

        ASSERT_LE(std::abs(mat.Elem(2,0) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(2,1) - 1./96), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - 5./96), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - 1./24), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - 1./16), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - 1./16), tol);
    }

    // check block(2,0)
    {
        auto mat = blockMassMatrix->GetBlock(2,0);

        ASSERT_LE(std::abs(mat.Elem(0,0) - 13./192), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(0,2) - 7./192), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - 1./192), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 1./32), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - 1./192), tol);

        ASSERT_LE(std::abs(mat.Elem(2,0) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(2,1) - 1./96), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - 5./96), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - 1./24), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - 1./16), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - 1./16), tol);

        ASSERT_LE(std::abs(mat.Elem(4,0) - 1./32), tol);
        ASSERT_LE(std::abs(mat.Elem(4,1) - 1./24), tol);
        ASSERT_LE(std::abs(mat.Elem(4,2) - 1./96), tol);
    }

    // check block(2,1)
    {
        auto mat = blockMassMatrix->GetBlock(2,1);

        ASSERT_LE(std::abs(mat.Elem(0,0) - 13./192), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - 1./192), tol);
        ASSERT_LE(std::abs(mat.Elem(0,2) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(0,3) - 1./32), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - 1./192), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 5./192), tol);
        ASSERT_LE(std::abs(mat.Elem(1,3) - 1./96), tol);

        ASSERT_LE(std::abs(mat.Elem(2,0) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - 1./24), tol);
        ASSERT_LE(std::abs(mat.Elem(2,3) - 1./48), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - 1./24), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - 1./48), tol);
        ASSERT_LE(std::abs(mat.Elem(3,3) - 1./12), tol);

        ASSERT_LE(std::abs(mat.Elem(4,0) - 1./32), tol);
        ASSERT_LE(std::abs(mat.Elem(4,1) - 1./32), tol);
        ASSERT_LE(std::abs(mat.Elem(4,3) - 1./48), tol);
    }

    // check upper-diagonal blocks
    {
        for (int n=1; n<3; n++)
            for (int m=0; m<n; m++)
            {
                auto matT = Transpose(blockMassMatrix->GetBlock(n,m));
                matT->Add(-1, blockMassMatrix->GetBlock(m,n));
                ASSERT_LE(matT->MaxNorm(), tol);
                delete matT;
            }
    }

    delete feColl;
}


/**
 * @brief Tests the assembly of sparse stiffness matrix
 */
TEST(SparseSpatialAssembly, stiffnessMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

    auto nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto stiffnessBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> stiffnessIntegator
            = std::make_shared<sparseHeat::SpatialStiffnessIntegrator>();

    stiffnessBilinearForm->addDomainIntegrator(stiffnessIntegator);
    stiffnessBilinearForm->assemble();

    auto blockStiffnessMatrix
            = stiffnessBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        BilinearForm stiffnessBf(fes1.get());
        stiffnessBf.AddDomainIntegrator(new DiffusionIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness00 = stiffnessBf.LoseMat();
        trueBlockStiffness00->Add(-1, blockStiffnessMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockStiffness00->MaxNorm(), tol);
    }
    {
        BilinearForm stiffnessBf(fes2.get());
        stiffnessBf.AddDomainIntegrator(new DiffusionIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness11 = stiffnessBf.LoseMat();
        trueBlockStiffness11->Add(-1, blockStiffnessMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockStiffness11->MaxNorm(), tol);
    }
    {
        BilinearForm stiffnessBf(fes3.get());
        stiffnessBf.AddDomainIntegrator(new DiffusionIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness22 = stiffnessBf.LoseMat();
        trueBlockStiffness22->Add(-1, blockStiffnessMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockStiffness22->MaxNorm(), tol);
    }

//    for(int i=1; i<blockStiffnessMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<i; j++) {
//            blockStiffnessMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    // check block(1,0)
    {
        auto mat = blockStiffnessMatrix->GetBlock(1,0);

        ASSERT_LE(std::abs(mat.Elem(0,0) - 1.), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (-1./2)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,2) - (-1./2)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,1) - (1./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - (-1./4)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (-1./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - (1./4)), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - (-1.)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - 1./2), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - 1./2), tol);
    }

    // check block(2,0)
    {
        auto mat = blockStiffnessMatrix->GetBlock(2,0);

        ASSERT_LE(std::abs(mat.Elem(0,0) - (3./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (-1./2)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,2) - (-1./4)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - (-1./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 1./4), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (-1./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - (1./4)), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - (-1.)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - 1./2), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - 1./2), tol);

        ASSERT_LE(std::abs(mat.Elem(4,0) - 1./2), tol);
        ASSERT_LE(std::abs(mat.Elem(4,2) - (-1./2)), tol);
    }

    // check block(2,1)
    {
        auto mat = blockStiffnessMatrix->GetBlock(2,1);

        ASSERT_LE(std::abs(mat.Elem(0,0) - 3./4), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (-1./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,3) - (-1./2)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - (-1./4)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - (1./4)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,2) - (1./2)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,3) - (-1./2)), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - (-1.)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - (-1./2)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - (-1./2)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,3) - (2.)), tol);

        ASSERT_LE(std::abs(mat.Elem(4,0) - 1./2), tol);
        ASSERT_LE(std::abs(mat.Elem(4,1) - 1./2), tol);
        ASSERT_LE(std::abs(mat.Elem(4,3) - (-1.)), tol);
    }

    // check upper-diagonal blocks
    {
        for (int n=1; n<3; n++)
            for (int m=0; m<n; m++)
            {
                auto matT
                        = Transpose(blockStiffnessMatrix->GetBlock(n,m));
                matT->Add(-1, blockStiffnessMatrix->GetBlock(m,n));
                ASSERT_LE(matT->MaxNorm(), tol);
                delete matT;
            }
    }

    delete feColl;
}


/**
 * @brief Tests the assembly of sparse VectorFE mass matrix
 */
TEST(SparseSpatialAssembly, vectorFEMassMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

    auto nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new RT_FECollection(0, dim);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto massBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> massIntegator
            = std::make_shared<sparseHeat::SpatialVectorFEMassIntegrator>();

    massBilinearForm->addDomainIntegrator(massIntegator);
    massBilinearForm->assemble();

    auto blockMassMatrix
            = massBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        BilinearForm massBf(fes1.get());
        massBf.AddDomainIntegrator(new VectorFEMassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass00 = massBf.LoseMat();
        trueBlockMass00->Add(-1, blockMassMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockMass00->MaxNorm(), tol);
    }
    {
        BilinearForm massBf(fes2.get());
        massBf.AddDomainIntegrator(new VectorFEMassIntegrator{});
        massBf.Assemble();
        auto trueBlockMass11 = massBf.LoseMat();
        trueBlockMass11->Add(-1, blockMassMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockMass11->MaxNorm(), tol);
    }
    {
        BilinearForm massBf(fes3.get());
        massBf.AddDomainIntegrator(new VectorFEMassIntegrator{});
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
 * @brief Tests the assembly of sparse VectorFE stiffness matrix
 */
TEST(SparseSpatialAssembly, vectorFEStiffnessMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

    auto nestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);

    int dim = mesh1->Dimension();
    FiniteElementCollection *feColl
            = new RT_FECollection(0, dim);
    auto fes1 = std::make_shared<FiniteElementSpace>(mesh1.get(), feColl);
    auto fes2 = std::make_shared<FiniteElementSpace>(mesh2.get(), feColl);
    auto fes3 = std::make_shared<FiniteElementSpace>(mesh3.get(), feColl);

    nestedFEHierarchy->addFESpace(fes1);
    nestedFEHierarchy->addFESpace(fes2);
    nestedFEHierarchy->addFESpace(fes3);

    auto stiffnessBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> stiffnessIntegator
            = std::make_shared<sparseHeat::SpatialVectorFEStiffnessIntegrator>();

    stiffnessBilinearForm->addDomainIntegrator(stiffnessIntegator);
    stiffnessBilinearForm->assemble();

    auto blockStiffnessMatrix = stiffnessBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        BilinearForm stiffnessBf(fes1.get());
        stiffnessBf.AddDomainIntegrator
                (new heat::VectorFEStiffnessIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness00 = stiffnessBf.LoseMat();
        trueBlockStiffness00->
                Add(-1, blockStiffnessMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockStiffness00->MaxNorm(), tol);
    }
    {
        BilinearForm stiffnessBf(fes2.get());
        stiffnessBf.AddDomainIntegrator
                (new heat::VectorFEStiffnessIntegrator{});
        stiffnessBf.Assemble();
        auto trueBlockStiffness11 = stiffnessBf.LoseMat();
        trueBlockStiffness11->
                Add(-1, blockStiffnessMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockStiffness11->MaxNorm(), tol);
    }
    {
        BilinearForm stiffnessBf(fes3.get());
        stiffnessBf.AddDomainIntegrator
                (new heat::VectorFEStiffnessIntegrator{});
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
 * @brief Tests the assembly of sparse VectorFE gradient matrix
 */
TEST(SparseSpatialAssembly, vectorFEGradientMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

    int dim = mesh1->Dimension();

    // Raviart-Thomas spaces
    FiniteElementCollection *feCollRT
            = new RT_FECollection(0, dim);
    auto fes1RT = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollRT);
    auto fes2RT = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollRT);
    auto fes3RT = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollRT);

    // H1 spaces
    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
    auto fes1H1 = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollH1);
    auto fes2H1 = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollH1);
    auto fes3H1 = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollH1);

    // hierarchy of test spaces
    auto testNestedFEHierarchy
            = std::make_shared<NestedFEHierarchy>(meshHierarchy);
    testNestedFEHierarchy->addFESpace(fes1RT);
    testNestedFEHierarchy->addFESpace(fes2RT);
    testNestedFEHierarchy->addFESpace(fes3RT);

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
            = std::make_shared<sparseHeat::SpatialVectorFEGradientIntegrator>();

    gradientBilinearForm->addDomainIntegrator(gradientIntegator);
    gradientBilinearForm->assemble();

    auto blockGradientMatrix
            = gradientBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        MixedBilinearForm gradientBf(fes1H1.get(),
                                     fes1RT.get());
        gradientBf.AddDomainIntegrator
                (new heat::VectorFEGradientIntegrator{});
        gradientBf.Assemble();
        auto trueBlockGradient00 = gradientBf.LoseMat();
        trueBlockGradient00->
                Add(-1, blockGradientMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockGradient00->MaxNorm(), tol);
    }
    {
        MixedBilinearForm gradientBf(fes2H1.get(),
                                     fes2RT.get());
        gradientBf.AddDomainIntegrator
                (new heat::VectorFEGradientIntegrator{});
        gradientBf.Assemble();
        auto trueBlockGradient11 = gradientBf.LoseMat();
        trueBlockGradient11->
                Add(-1, blockGradientMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockGradient11->MaxNorm(), tol);
    }
    {
        MixedBilinearForm gradientBf(fes3H1.get(),
                                     fes3RT.get());
        gradientBf.AddDomainIntegrator
                (new heat::VectorFEGradientIntegrator{});
        gradientBf.Assemble();
        auto trueBlockGradient22 = gradientBf.LoseMat();
        trueBlockGradient22->
                Add(-1, blockGradientMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockGradient22->MaxNorm(), tol);
    }

    for(int i=1; i<blockGradientMatrix->NumRowBlocks(); i++) {
        for (int j=0; j<i; j++) {
            blockGradientMatrix->GetBlock(i,j).Print();
            std::cout << "\n";
        }
    }

    delete feCollRT;
    delete feCollH1;
}


/**
 * @brief Tests the assembly of sparse VectorFE divergence matrix
 */
TEST(SparseSpatialAssembly, vectorFEDivergenceMatrix)
{
    std::string input_dir
            = "../tests/input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

    int dim = mesh1->Dimension();

    // Raviart-Thomas spaces
    FiniteElementCollection *feCollRT
            = new RT_FECollection(0, dim);
    auto fes1RT = std::make_shared<FiniteElementSpace>(mesh1.get(),
                                                       feCollRT);
    auto fes2RT = std::make_shared<FiniteElementSpace>(mesh2.get(),
                                                       feCollRT);
    auto fes3RT = std::make_shared<FiniteElementSpace>(mesh3.get(),
                                                       feCollRT);

    // H1 spaces
    FiniteElementCollection *feCollH1
            = new H1_FECollection(1, dim, BasisType::GaussLobatto);
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
    trialNestedFEHierarchy->addFESpace(fes1RT);
    trialNestedFEHierarchy->addFESpace(fes2RT);
    trialNestedFEHierarchy->addFESpace(fes3RT);

    // assemble divergence mixed bilinear form
    auto divergenceBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (trialNestedFEHierarchy, testNestedFEHierarchy);

    std::shared_ptr<BlockMixedBilinearFormIntegrator> divergenceIntegator
            = std::make_shared<sparseHeat::SpatialVectorFEDivergenceIntegrator>();

    divergenceBilinearForm->addDomainIntegrator(divergenceIntegator);
    divergenceBilinearForm->assemble();

    auto blockDivergenceMatrix
            = divergenceBilinearForm->getBlockMatrix();

    // check diagonal blocks
    double tol = 1E-8;
    {
        MixedBilinearForm divergenceBf(fes1RT.get(),
                                       fes1H1.get());
        divergenceBf.AddDomainIntegrator
                (new VectorFEDivergenceIntegrator{});
        divergenceBf.Assemble();
        auto trueBlockDivergence00 = divergenceBf.LoseMat();
        trueBlockDivergence00->
                Add(-1, blockDivergenceMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockDivergence00->MaxNorm(), tol);
    }
    {
        MixedBilinearForm divergenceBf(fes2RT.get(),
                                       fes2H1.get());
        divergenceBf.AddDomainIntegrator
                (new VectorFEDivergenceIntegrator{});
        divergenceBf.Assemble();
        auto trueBlockDivergence11 = divergenceBf.LoseMat();
        trueBlockDivergence11->
                Add(-1, blockDivergenceMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockDivergence11->MaxNorm(), tol);
    }
    {
        MixedBilinearForm divergenceBf(fes3RT.get(),
                                       fes3H1.get());
        divergenceBf.AddDomainIntegrator
                (new VectorFEDivergenceIntegrator{});
        divergenceBf.Assemble();
        auto trueBlockDivergence22 = divergenceBf.LoseMat();
        trueBlockDivergence22->
                Add(-1, blockDivergenceMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockDivergence22->MaxNorm(), tol);
    }

    for(int i=1; i<blockDivergenceMatrix->NumRowBlocks(); i++) {
        for (int j=0; j<i; j++) {
            blockDivergenceMatrix->GetBlock(i,j).Print();
            std::cout << "\n";
        }
    }

    delete feCollRT;
    delete feCollH1;
}
