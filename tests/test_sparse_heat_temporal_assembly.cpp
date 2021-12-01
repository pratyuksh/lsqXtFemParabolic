#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../src/core/config.hpp"
#include "../src/heat/assembly.hpp"
#include "../src/sparse_heat/temporal_assembly.hpp"


/**
 * @brief Tests the assembly of sparse mass matrix
 */
TEST(SparseTemporalAssembly, massMatrix)
{
    double endTime = 1;

    int minLevel = 1;
    int maxLevel = 4;
    int numLevels = 4;

    auto massAssembler
            = std::make_unique<sparseHeat::TemporalMassMatrixAssembler>
            (endTime, minLevel, maxLevel);

    auto blockSizes = massAssembler->getBlockSizes();
    auto blockOffsets = massAssembler->getBlockOffsets();

    ASSERT_EQ(blockSizes.Size(), numLevels);
    ASSERT_EQ(blockOffsets.Size(), numLevels+1);

    Array<int> trueBlockSizes(numLevels);
    trueBlockSizes[0] = 3;
    trueBlockSizes[1] = 2;
    trueBlockSizes[2] = 4;
    trueBlockSizes[3] = 8;
    ASSERT_EQ(blockSizes, trueBlockSizes);

    Array<int> trueBlockOffsets(numLevels+1);
    trueBlockOffsets[0] = 0;
    trueBlockOffsets[1] = 3;
    trueBlockOffsets[2] = 5;
    trueBlockOffsets[3] = 9;
    trueBlockOffsets[4] = 17;
    ASSERT_EQ(blockOffsets, trueBlockOffsets);

    massAssembler->assemble();
    auto blockMassMatrix = massAssembler->getBlockMatrix();

//    for(int i=0; i<blockMassMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<=i; j++) {
//            blockMassMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    // check diagonal blocks
    double tol = 1E-8;
    {
        auto mat = blockMassMatrix->GetBlock(0, 0);
        double h = endTime / pow(2, minLevel);

        ASSERT_LE(std::abs(mat.Elem(0,0) - h/3), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - h/6), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - h/6), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 2*h/3), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - h/6), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - h/6), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - h/3), tol);
    }
    for (int m = 1; m < numLevels; m++)
    {
        int level = minLevel + m;
        auto mat = blockMassMatrix->GetBlock(m, m);
        double h = endTime / std::pow(2, level);

        for (int i=0; i<blockSizes[m]; i++) {
            ASSERT_LE(std::abs(mat.Elem(i,i) - 2*h/3), tol);
        }
    }

    // check lower-diagonal blocks
    { // block(1, 0)
        auto mat = blockMassMatrix->GetBlock(1, 0);
        double h0 = endTime / pow(2, minLevel);
        double h1 = h0/2;
        double h1SquaredDivH0 = h1 * h1 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - h1SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - h1SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(1,1) - h1SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - h1SquaredDivH0), tol);
    }
    { // block(2, 1)
        auto mat = blockMassMatrix->GetBlock(2, 1);
        double h1 = endTime / pow(2, minLevel+1);
        double h2 = h1/2;
        double h2SquaredDivH1 = h2 * h2 / h1;

        ASSERT_LE(std::abs(mat.Elem(0,0) - h2SquaredDivH1), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - h2SquaredDivH1), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - h2SquaredDivH1), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - h2SquaredDivH1), tol);
    }
    { // block(2, 0)
        auto mat = blockMassMatrix->GetBlock(2, 0);
        double h0 = endTime / pow(2, minLevel);
        double h2 = h0/4;
        double h2SquaredDivH0 = h2 * h2 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - 3*h2SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - h2SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - h2SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 3*h2SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - 3*h2SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - h2SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(3,1) - h2SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - 3*h2SquaredDivH0), tol);
    }
    { // block(3, 2)
        auto mat = blockMassMatrix->GetBlock(3, 2);
        double h2 = endTime / pow(2, minLevel+2);
        double h3 = h2/2;
        double h3SquaredDivH2 = h3 * h3 / h2;

        ASSERT_LE(std::abs(mat.Elem(0,0) - h3SquaredDivH2), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - h3SquaredDivH2), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - h3SquaredDivH2), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - h3SquaredDivH2), tol);

        ASSERT_LE(std::abs(mat.Elem(4,2) - h3SquaredDivH2), tol);
        ASSERT_LE(std::abs(mat.Elem(5,2) - h3SquaredDivH2), tol);

        ASSERT_LE(std::abs(mat.Elem(6,3) - h3SquaredDivH2), tol);
        ASSERT_LE(std::abs(mat.Elem(7,3) - h3SquaredDivH2), tol);
    }
    { // block(3, 1)
        auto mat = blockMassMatrix->GetBlock(3, 1);
        double h1 = endTime / pow(2, minLevel+1);
        double h3 = h1/4;
        double h3SquaredDivH1 = h3 * h3 / h1;

        ASSERT_LE(std::abs(mat.Elem(0,0) - h3SquaredDivH1), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - 3*h3SquaredDivH1), tol);

        ASSERT_LE(std::abs(mat.Elem(2,0) - 3*h3SquaredDivH1), tol);
        ASSERT_LE(std::abs(mat.Elem(3,0) - h3SquaredDivH1), tol);

        ASSERT_LE(std::abs(mat.Elem(4,1) - h3SquaredDivH1), tol);
        ASSERT_LE(std::abs(mat.Elem(5,1) - 3*h3SquaredDivH1), tol);

        ASSERT_LE(std::abs(mat.Elem(6,1) - 3*h3SquaredDivH1), tol);
        ASSERT_LE(std::abs(mat.Elem(7,1) - h3SquaredDivH1), tol);
    }
    { // block(3, 0)
        auto mat = blockMassMatrix->GetBlock(3, 0);
        double h0 = endTime / pow(2, minLevel);
        double h3 = h0/8;
        double h3SquaredDivH0 = h3 * h3 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - 7*h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - 5*h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - 3*h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(2,0) - 3*h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(2,1) - 5*h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - 7*h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(4,1) - 7*h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(4,2) - h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(5,1) - 5*h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(5,2) - 3*h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(6,1) - 3*h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(6,2) - 5*h3SquaredDivH0), tol);

        ASSERT_LE(std::abs(mat.Elem(7,1) - h3SquaredDivH0), tol);
        ASSERT_LE(std::abs(mat.Elem(7,2) - 7*h3SquaredDivH0), tol);
    }

    // check upper-diagonal blocks
    {
        for (int n=1; n < numLevels; n++)
            for (int m=0; m < n; m++)
            {
                auto matT = Transpose
                        (blockMassMatrix->GetBlock(n,m));
                matT->Add(-1, blockMassMatrix->GetBlock(m,n));
                ASSERT_LE(matT->MaxNorm(), tol);
                delete matT;
            }
    }
}


/**
 * @brief Tests the assembly of sparse stiffness matrix
 */
TEST(SparseTemporalAssembly, stiffnessMatrix)
{
    double endTime = 1;

    int minLevel = 1;
    int maxLevel = 4;
    int numLevels = 4;

    auto stiffnessAssembler
            = std::make_unique<sparseHeat::TemporalStiffnessMatrixAssembler>
            (endTime, minLevel, maxLevel);

    auto blockSizes = stiffnessAssembler->getBlockSizes();
    auto blockOffsets = stiffnessAssembler->getBlockOffsets();

    ASSERT_EQ(blockSizes.Size(), numLevels);
    ASSERT_EQ(blockOffsets.Size(), numLevels+1);

    Array<int> trueBlockSizes(numLevels);
    trueBlockSizes[0] = 3;
    trueBlockSizes[1] = 2;
    trueBlockSizes[2] = 4;
    trueBlockSizes[3] = 8;
    ASSERT_EQ(blockSizes, trueBlockSizes);

    Array<int> trueBlockOffsets(numLevels+1);
    trueBlockOffsets[0] = 0;
    trueBlockOffsets[1] = 3;
    trueBlockOffsets[2] = 5;
    trueBlockOffsets[3] = 9;
    trueBlockOffsets[4] = 17;
    ASSERT_EQ(blockOffsets, trueBlockOffsets);

    stiffnessAssembler->assemble();
    auto blockStiffnessMatrix = stiffnessAssembler->getBlockMatrix();

//    for(int i=0; i<blockStiffnessMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<=i; j++) {
//            blockStiffnessMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    // check diagonal blocks
    double tol = 1E-8;
    {
        auto mat = blockStiffnessMatrix->GetBlock(0, 0);
        double h = endTime / pow(2, minLevel);

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+1/h)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (-1/h)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - (-1/h)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - (+2/h)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - (-1/h)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (-1/h)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - (+1/h)), tol);
    }
    for (int m = 1; m < numLevels; m++)
    {
        int level = minLevel + m;
        auto mat = blockStiffnessMatrix->GetBlock(m, m);
        double h = endTime / std::pow(2, level);

        for (int i=0; i<blockSizes[m]; i++) {
            ASSERT_LE(std::abs(mat.Elem(i,i) - 2/h), tol);
        }
    }

    // check off-diagonal blocks, they are zero
    {
        for (int m=1; m < numLevels; m++)
            for (int n=0; n < m; n++)
            {
                auto mat1 = blockStiffnessMatrix->GetBlock(m,n);
                auto mat2 = blockStiffnessMatrix->GetBlock(n,m);
                ASSERT_LE(mat1.MaxNorm(), tol);
                ASSERT_LE(mat2.MaxNorm(), tol);
            }
    }
}


/**
 * @brief Tests the assembly of sparse gradient matrix
 */
TEST(SparseTemporalAssembly, gradientMatrix)
{
    double endTime = 1;

    int minLevel = 1;
    int maxLevel = 4;
    int numLevels = 4;

    auto gradientAssembler
            = std::make_unique<sparseHeat::TemporalGradientMatrixAssembler>
            (endTime, minLevel, maxLevel);

    auto blockSizes = gradientAssembler->getBlockSizes();
    auto blockOffsets = gradientAssembler->getBlockOffsets();

    ASSERT_EQ(blockSizes.Size(), numLevels);
    ASSERT_EQ(blockOffsets.Size(), numLevels+1);

    Array<int> trueBlockSizes(numLevels);
    trueBlockSizes[0] = 3;
    trueBlockSizes[1] = 2;
    trueBlockSizes[2] = 4;
    trueBlockSizes[3] = 8;
    ASSERT_EQ(blockSizes, trueBlockSizes);

    Array<int> trueBlockOffsets(numLevels+1);
    trueBlockOffsets[0] = 0;
    trueBlockOffsets[1] = 3;
    trueBlockOffsets[2] = 5;
    trueBlockOffsets[3] = 9;
    trueBlockOffsets[4] = 17;
    ASSERT_EQ(blockOffsets, trueBlockOffsets);

    gradientAssembler->assemble();
    auto blockGradientMatrix = gradientAssembler->getBlockMatrix();

//    for(int i=0; i<blockGradientMatrix->NumRowBlocks(); i++) {
//        for (int j=0; j<=i; j++) {
//            blockGradientMatrix->GetBlock(i,j).Print();
//            std::cout << "\n";
//        }
//    }

    // check diagonal blocks
    double tol = 1E-8;
    {
        auto mat = blockGradientMatrix->GetBlock(0, 0);
        double half = 1./2;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-half)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (+half)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - (-half)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - (+half)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (-half)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - (+half)), tol);
    }
    for (int m = 1; m < numLevels; m++)
    {
        auto mat = blockGradientMatrix->GetBlock(m, m);
        for (int i=0; i<blockSizes[m]; i++) {
            ASSERT_LE(mat.MaxNorm(), tol);
        }
    }

    // check lower-diagonal blocks
    { // block(1, 0)
        auto mat = blockGradientMatrix->GetBlock(1, 0);
        double h0 = endTime / pow(2, minLevel);
        double h1 = h0/2;
        double h1DivH0 = h1 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-h1DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h1DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,1) - (-h1DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - (+h1DivH0)), tol);
    }
    { // block(2, 1)
        auto mat = blockGradientMatrix->GetBlock(2, 1);
        double h1 = endTime / pow(2, minLevel+1);
        double h2 = h1/2;
        double h2DivH1 = h2 / h1;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+h2DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h2DivH1)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (+h2DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - (-h2DivH1)), tol);
    }
    { // block(2, 0)
        auto mat = blockGradientMatrix->GetBlock(2, 0);
        double h0 = endTime / pow(2, minLevel);
        double h2 = h0/4;
        double h2DivH0 = h2 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h2DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - (+h2DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (-h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - (+h2DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(3,1) - (-h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,2) - (+h2DivH0)), tol);
    }
    { // block(3, 2)
        auto mat = blockGradientMatrix->GetBlock(3, 2);
        double h2 = endTime / pow(2, minLevel+2);
        double h3 = h2/2;
        double h3DivH2 = h3 / h2;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h3DivH2)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,1) - (+h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - (-h3DivH2)), tol);

        ASSERT_LE(std::abs(mat.Elem(4,2) - (+h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(5,2) - (-h3DivH2)), tol);

        ASSERT_LE(std::abs(mat.Elem(6,3) - (+h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(7,3) - (-h3DivH2)), tol);
    }
    { // block(3, 1)
        auto mat = blockGradientMatrix->GetBlock(3, 1);
        double h1 = endTime / pow(2, minLevel+1);
        double h3 = h1/4;
        double h3DivH1 = h3 / h1;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - (+h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,0) - (-h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,0) - (-h3DivH1)), tol);

        ASSERT_LE(std::abs(mat.Elem(4,1) - (+h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(5,1) - (+h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(6,1) - (-h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(7,1) - (-h3DivH1)), tol);
    }
    { // block(3, 0)
        auto mat = blockGradientMatrix->GetBlock(3, 0);
        double h0 = endTime / pow(2, minLevel);
        double h3 = h0/8;
        double h3DivH0 = h3 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,0) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,1) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(3,0) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,1) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(4,1) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(4,2) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(5,1) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(5,2) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(6,1) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(6,2) - (+h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(7,1) - (-h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(7,2) - (+h3DivH0)), tol);
    }

    // check pper-diagonal blocks
    { // block(0, 1)
        auto mat = blockGradientMatrix->GetBlock(0, 1);
        double h0 = endTime / pow(2, minLevel);
        double h1 = h0/2;
        double h1DivH0 = h1 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+h1DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h1DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,1) - (+h1DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,1) - (-h1DivH0)), tol);
    }
    { // block(1, 2)
        auto mat = blockGradientMatrix->GetBlock(1, 2);
        double h1 = endTime / pow(2, minLevel+1);
        double h2 = h1/2;
        double h2DivH1 = h2 / h1;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-h2DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h2DivH1)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,2) - (-h2DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,3) - (+h2DivH1)), tol);
    }
    { // block(0, 2)
        auto mat = blockGradientMatrix->GetBlock(0, 2);
        double h0 = endTime / pow(2, minLevel);
        double h2 = h0/4;
        double h2DivH0 = h2 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h2DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - (-h2DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,2) - (+h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,2) - (-h2DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,3) - (+h2DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,3) - (-h2DivH0)), tol);
    }
    { // block(2, 3)
        auto mat = blockGradientMatrix->GetBlock(2, 3);
        double h2 = endTime / pow(2, minLevel+2);
        double h3 = h2/2;
        double h3DivH2 = h3 / h2;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h3DivH2)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,2) - (-h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,3) - (+h3DivH2)), tol);

        ASSERT_LE(std::abs(mat.Elem(2,4) - (-h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,5) - (+h3DivH2)), tol);

        ASSERT_LE(std::abs(mat.Elem(3,6) - (-h3DivH2)), tol);
        ASSERT_LE(std::abs(mat.Elem(3,7) - (+h3DivH2)), tol);
    }
    { // block(1, 3)
        auto mat = blockGradientMatrix->GetBlock(1, 3);
        double h1 = endTime / pow(2, minLevel+1);
        double h3 = h1/4;
        double h3DivH1 = h3 / h1;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (-h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,1) - (-h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,2) - (+h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(0,3) - (+h3DivH1)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,4) - (-h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,5) - (-h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,6) - (+h3DivH1)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,7) - (+h3DivH1)), tol);
    }
    { // block(0, 3)
        auto mat = blockGradientMatrix->GetBlock(0, 3);
        double h0 = endTime / pow(2, minLevel);
        double h3 = h0/8;
        double h3DivH0 = h3 / h0;

        ASSERT_LE(std::abs(mat.Elem(0,0) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,0) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(0,1) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,1) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(0,2) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,2) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(0,3) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(1,3) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,4) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,4) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,5) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,5) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,6) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,6) - (-h3DivH0)), tol);

        ASSERT_LE(std::abs(mat.Elem(1,7) - (+h3DivH0)), tol);
        ASSERT_LE(std::abs(mat.Elem(2,7) - (-h3DivH0)), tol);
    }
}
