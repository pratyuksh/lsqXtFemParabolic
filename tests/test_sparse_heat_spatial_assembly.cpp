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

// End of file
