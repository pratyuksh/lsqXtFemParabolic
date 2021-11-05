#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>
#include <chrono>

#include "../src/core/config.hpp"
#include "../src/sparse_heat/assembly.hpp"


using namespace mymfem;

/**
 * @brief Tests the assembly of sparse mass matrix
 */
TEST(SparseAssembly, massMatrix)
{
    std::string input_dir
            = "../input/sparse_assembly/";

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
    meshHierarchy->buildHierarchicalTranformations();

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

    std::unique_ptr<BlockBilinearForm> massBilinearForm
            = std::make_unique<BlockBilinearForm>(nestedFEHierarchy);

    std::shared_ptr<BlockBilinearFormIntegrator> massIntegator
            = std::make_shared<sparseHeat::MassIntegrator>();

    massBilinearForm->addDomainIntegrator(massIntegator);
    massBilinearForm->assemble();

    auto blockMassMatrix = massBilinearForm->getBlockMatrix();

    // check diagonal blocks
    {
        BilinearForm massBf(fes1.get());
        massBf.AddDomainIntegrator(new MassIntegrator);
        massBf.Assemble();
        auto trueBlockMass00 = massBf.LoseMat();
        trueBlockMass00->Add(-1, blockMassMatrix->GetBlock(0,0));
        ASSERT_LE(trueBlockMass00->MaxNorm(), 1E-8);
    }
    {
        BilinearForm massBf(fes2.get());
        massBf.AddDomainIntegrator(new MassIntegrator);
        massBf.Assemble();
        auto trueBlockMass11 = massBf.LoseMat();
        trueBlockMass11->Add(-1, blockMassMatrix->GetBlock(1,1));
        ASSERT_LE(trueBlockMass11->MaxNorm(), 1E-8);
    }
    {
        BilinearForm massBf(fes3.get());
        massBf.AddDomainIntegrator(new MassIntegrator);
        massBf.Assemble();
        auto trueBlockMass22 = massBf.LoseMat();
        trueBlockMass22->Add(-1, blockMassMatrix->GetBlock(2,2));
        ASSERT_LE(trueBlockMass22->MaxNorm(), 1E-8);
    }

    /*for(int i=1; i<blockMassMatrix->NumRowBlocks(); i++) {
        for (int j=0; j<i; j++) {
            blockMassMatrix->GetBlock(i,j).Print();
            std::cout << "\n";
        }
    }

    for(int i=0; i<blockMassMatrix->NumRowBlocks(); i++) {
        for (int j=i+1; j<blockMassMatrix->NumColBlocks(); j++) {
            blockMassMatrix->GetBlock(i,j).Print();
            std::cout << "\n";
        }
    }*/
}
