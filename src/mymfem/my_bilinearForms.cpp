#include "my_bilinearForms.hpp"
#include "utilities.hpp"

using namespace mfem;


//! Constructor with nested finite element hierarchy
//! defines the block offsets of storage matrix
mymfem::BlockBilinearForm
:: BlockBilinearForm
(const std::shared_ptr<NestedFEHierarchy>& nestedFEHierarchy)
    : m_nestedFEHierarchy (nestedFEHierarchy)
{
    m_numLevels = m_nestedFEHierarchy->getNumLevels();
    auto feSpaces = m_nestedFEHierarchy->getFESpaces();

    // define block offsets from FE spaces
    m_blockOffsets.SetSize(m_numLevels+1);
    m_blockOffsets[0] = 0;
    for (int i=0; i<m_numLevels; i++) {
        m_blockOffsets[i+1] = feSpaces[i]->GetTrueVSize();
    }
    m_blockOffsets.PartialSum();
}

//! Assembles the domain block bilinear form integrators
void mymfem::BlockBilinearForm
:: assemble (bool symmetric, int skip_zeros)
{
    // if not initialized, allocate block matrix
    if (!m_blockMatrix) {
        this->allocateBlockMatrix();
    }

    if (mydbfi.size())
    {
        auto feSpaces = m_nestedFEHierarchy->getFESpaces();
        auto hierMeshTrans
                = m_nestedFEHierarchy->getHierarchicalMeshTransformations();

        DenseMatrix elmat;        

        // Diagonal blocks
        Array<int> vdofs;
        ElementTransformation *elTrans = nullptr;

        for (int n=0; n<m_numLevels; n++)
        {
            auto& mat = m_blockMatrix->GetBlock(n,n);
            std::shared_ptr<FiniteElementSpace> fes(feSpaces[n]);

            for (int i = 0; i < fes -> GetNE(); i++)
            {
                fes->GetElementVDofs(i, vdofs);

                const FiniteElement &fe = *(fes->GetFE(i));
                elTrans = fes->GetElementTransformation(i);
                for (int k = 0; k < static_cast<int>(mydbfi.size()); k++)
                {
                    mydbfi[k]->assembleElementMatrix
                            (fe, *elTrans, elmat);
                    mat.AddSubMatrix(vdofs, vdofs, elmat, skip_zeros);
                }
            }            
        }

        // Lower-triangular blocks
        Array<int> trialVdofsCoarse, testVdofsFine;
        ElementTransformation *trialElTransCoarse = nullptr;
        ElementTransformation *testElTransFine = nullptr;

        for (int m=1; m<m_numLevels; m++)
        {
            std::shared_ptr<FiniteElementSpace> fesFine(feSpaces[m]);
            auto hierMultiLevelMeshTrans = hierMeshTrans[m-1];

            for (int n=m-1; n>=0; n--)
            {
                std::shared_ptr<FiniteElementSpace> fesCoarse(feSpaces[n]);
                // generate hierarchical mesh transformations
                // across non-successive mesh levels
                if (n < m-1) {
                    hierMultiLevelMeshTrans
                            = buildMultiLevelMeshHierarchy
                            (hierMeshTrans[n], hierMultiLevelMeshTrans);
                }

                auto& mat = m_blockMatrix->GetBlock(m, n);
                for (int i=0; i < fesCoarse->GetNE(); i++)
                {
                    fesCoarse->GetElementVDofs(i, trialVdofsCoarse);
                    const FiniteElement &trialFECoarse
                            = *(fesCoarse->GetFE(i));
                    trialElTransCoarse
                            = fesCoarse->GetElementTransformation(i);

                    auto childElems = (*hierMultiLevelMeshTrans)(i);
                    for (auto it=childElems.begin();
                         it!= childElems.end(); it++)
                    {
                        int j = *it;
//                        std::cout << n << "\t" << m << "\t"
//                                  << i << "\t" << j << std::endl;

                        fesFine->GetElementVDofs(j, testVdofsFine);
                        const FiniteElement &testFEFine
                                = *(fesFine->GetFE(j));
                        testElTransFine
                                = fesFine->GetElementTransformation(j);

                        // call the integrators here
                        for (int k = 0; k < static_cast<int>(mydbfi.size()); k++)
                        {
                            mydbfi[k]->assembleElementMatrix
                                    (testFEFine, *testElTransFine,
                                     trialFECoarse, *trialElTransCoarse,
                                     elmat);
                            mat.AddSubMatrix(testVdofsFine,
                                             trialVdofsCoarse,
                                             elmat, skip_zeros);
                        }
                    }
                }
//                std::cout << "\n\n";
            }
        }
        m_blockMatrix->Finalize(skip_zeros, false);

        // upper-triangle blocks
        if (symmetric) {
            for (int n=1; n<m_numLevels; n++)
                for (int m=0; m<n; m++)
                {
                    auto matT = Transpose(m_blockMatrix->GetBlock(n,m));
                    m_blockMatrix->GetBlock(m,n)
                            = *matT;
                    delete matT;
                }
        }
        else {
            //TODO: compute upper-triangular blocks for non-symmetric case
        }
    }
}

void mymfem::BlockBilinearForm
:: allocateBlockMatrix()
{
    auto feSpaces = m_nestedFEHierarchy->getFESpaces();

    m_blockMatrix = std::make_shared<BlockMatrix>(m_blockOffsets);
    for(int i=0; i<m_numLevels; i++) {
        for(int j=0; j<m_numLevels; j++) {
            SparseMatrix *mat
                    = new SparseMatrix(feSpaces[i]->GetTrueVSize(),
                                       feSpaces[j]->GetTrueVSize());
            m_blockMatrix->SetBlock(i, j, mat);
        }
    }

    // owns the allocated memory
    m_blockMatrix->owns_blocks = true;
}


//! Constructor with nested finite element hierarchy
//! defines the block offsets of storage matrix
mymfem::BlockMixedBilinearForm
:: BlockMixedBilinearForm
(const std::shared_ptr<NestedFEHierarchy>& trialNestedFEHierarchy,
 const std::shared_ptr<NestedFEHierarchy>& testNestedFEHierarchy)
    : m_trialNestedFEHierarchy (trialNestedFEHierarchy),
      m_testNestedFEHierarchy (testNestedFEHierarchy)
{
    m_numLevels = m_testNestedFEHierarchy->getNumLevels();
    assert(m_numLevels == m_trialNestedFEHierarchy->getNumLevels());

    auto testFeSpaces = m_testNestedFEHierarchy->getFESpaces();
    auto trialFeSpaces = m_trialNestedFEHierarchy->getFESpaces();

    // define row block offsets from test FE spaces
    m_rowBlockOffsets.SetSize(m_numLevels+1);
    m_rowBlockOffsets[0] = 0;
    for (int i=0; i<m_numLevels; i++) {
        m_rowBlockOffsets[i+1] = testFeSpaces[i]->GetTrueVSize();
    }
    m_rowBlockOffsets.PartialSum();

    // define column block offsets from trial FE spaces
    m_colBlockOffsets.SetSize(m_numLevels+1);
    m_colBlockOffsets[0] = 0;
    for (int i=0; i<m_numLevels; i++) {
        m_colBlockOffsets[i+1] = trialFeSpaces[i]->GetTrueVSize();
    }
    m_colBlockOffsets.PartialSum();
}

//! Assembles the domain block bilinear form integrators
void mymfem::BlockMixedBilinearForm
:: assemble (int skip_zeros)
{
    // if not initialized, allocate block matrix
    if (!m_blockMatrix) {
        this->allocateBlockMatrix();
    }

    if (mydbfi.size())
    {
        // The underlying mesh hierarchy is the same for both test and trial spaces
        auto hierMeshTrans
                = m_testNestedFEHierarchy->getHierarchicalMeshTransformations();

        auto testFeSpaces = m_testNestedFEHierarchy->getFESpaces();
        auto trialFeSpaces = m_trialNestedFEHierarchy->getFESpaces();

        // Diagonal blocks
        DenseMatrix elmat;
        Array<int> testVdofs, trialVdofs;
        ElementTransformation *elTrans = nullptr;

        for (int n=0; n<m_numLevels; n++)
        {
            auto& mat = m_blockMatrix->GetBlock(n,n);
            std::shared_ptr<FiniteElementSpace> testFes(testFeSpaces[n]);
            std::shared_ptr<FiniteElementSpace> trialFes(trialFeSpaces[n]);

            for (int i = 0; i < testFes -> GetNE(); i++)
            {
                testFes->GetElementVDofs(i, testVdofs);
                trialFes->GetElementVDofs(i, trialVdofs);

                const FiniteElement &testFe = *(testFes->GetFE(i));
                const FiniteElement &trialFe = *(trialFes->GetFE(i));
                elTrans = testFes->GetElementTransformation(i);
                for (int k = 0; k < static_cast<int>(mydbfi.size()); k++)
                {
                    mydbfi[k]->assembleElementMatrix
                            (testFe, trialFe, *elTrans, elmat);
                    mat.AddSubMatrix(testVdofs, trialVdofs, elmat, skip_zeros);
                }
            }
        }

        // Lower-triangular and upper-triangular blocks
        DenseMatrix elmat1, elmat2;

        // auxiliary variables for lower-triangular blocks
        Array<int> trialVdofsCoarse, testVdofsFine;
        ElementTransformation *trialElTransCoarse = nullptr;
        ElementTransformation *testElTransFine = nullptr;
        std::shared_ptr<FiniteElementSpace> trialFesCoarse, testFesFine;

        // auxiliary variables for upper-triangular blocks
        Array<int> trialVdofsFine, testVdofsCoarse;
        ElementTransformation *trialElTransFine = nullptr;
        ElementTransformation *testElTransCoarse = nullptr;
        std::shared_ptr<FiniteElementSpace> trialFesFine, testFesCoarse;

        for (int m=1; m<m_numLevels; m++)
        {
            auto hierMultiLevelMeshTrans = hierMeshTrans[m-1];

            testFesFine = testFeSpaces[m];
            trialFesFine = trialFeSpaces[m];

            for (int n=m-1; n>=0; n--)
            {
                testFesCoarse = testFeSpaces[n];
                trialFesCoarse = trialFeSpaces[n];

                // generate hierarchical mesh transformations
                // across non-successive mesh levels
                if (n < m-1) {
                    hierMultiLevelMeshTrans
                            = buildMultiLevelMeshHierarchy
                            (hierMeshTrans[n], hierMultiLevelMeshTrans);
                }

                auto& mat1 = m_blockMatrix->GetBlock(m, n);
                auto& mat2 = m_blockMatrix->GetBlock(n, m);
                for (int i=0; i < trialFesCoarse->GetNE(); i++)
                {
                    // set variables for computing the matrices in lower-triangular blocks
                    trialFesCoarse->GetElementVDofs(i, trialVdofsCoarse);
                    const FiniteElement &trialFECoarse
                            = *(trialFesCoarse->GetFE(i));
                    trialElTransCoarse
                            = trialFesCoarse->GetElementTransformation(i);

                    // set variables for computing the matrices in upper-triangular blocks
                    testFesCoarse->GetElementVDofs(i, testVdofsCoarse);
                    const FiniteElement &testFECoarse
                            = *(testFesCoarse->GetFE(i));
                    testElTransCoarse
                            = testFesCoarse->GetElementTransformation(i);

                    auto childElems = (*hierMultiLevelMeshTrans)(i);
                    for (auto it=childElems.begin();
                         it!= childElems.end(); it++)
                    {
                        int j = *it;
//                        std::cout << n << "\t" << m << "\t"
//                                  << i << "\t" << j << std::endl;

                        // set variables for computing the matrices in the lower-triangular blocks
                        testFesFine->GetElementVDofs(j, testVdofsFine);
                        const FiniteElement &testFEFine
                                = *(testFesFine->GetFE(j));
                        testElTransFine
                                = testFesFine->GetElementTransformation(j);

                        // set variables for computing the matrices in the upper-triangular blocks
                        trialFesFine->GetElementVDofs(j, trialVdofsFine);
                        const FiniteElement &trialFEFine
                                = *(trialFesFine->GetFE(j));
                        trialElTransFine
                                = trialFesFine->GetElementTransformation(j);

                        // call the integrators here
                        for (int k = 0; k < static_cast<int>(mydbfi.size()); k++)
                        {
                            // compute matrices in the lower-triangular blocks
                            mydbfi[k]->assembleElementMatrix
                                    (testFEFine, *testElTransFine,
                                     trialFECoarse, *trialElTransCoarse,
                                     elmat1);
                            mat1.AddSubMatrix(testVdofsFine,
                                              trialVdofsCoarse,
                                              elmat1, skip_zeros);

                            // compute matrices in the uppwer-triangular blocks
                            mydbfi[k]->assembleElementMatrix2
                                    (testFECoarse, *testElTransCoarse,
                                     trialFEFine, *trialElTransFine,
                                     elmat2);
                            mat2.AddSubMatrix(testVdofsCoarse,
                                              trialVdofsFine,
                                              elmat2, skip_zeros);
                        }
                    }
                }
//                std::cout << "\n\n";
            }
        }
        m_blockMatrix->Finalize(skip_zeros, false);
    }
}

void mymfem::BlockMixedBilinearForm
:: allocateBlockMatrix()
{
    auto testFeSpaces = m_testNestedFEHierarchy->getFESpaces();
    auto trialFeSpaces = m_trialNestedFEHierarchy->getFESpaces();

    m_blockMatrix = std::make_shared<BlockMatrix>(m_rowBlockOffsets,
                                                  m_colBlockOffsets);
    for(int i=0; i<m_numLevels; i++) {
        for(int j=0; j<m_numLevels; j++) {
            SparseMatrix *mat
                    = new SparseMatrix(testFeSpaces[i]->GetTrueVSize(),
                                       trialFeSpaces[j]->GetTrueVSize());
            m_blockMatrix->SetBlock(i, j, mat);
        }
    }

    // owns the allocated memory
    m_blockMatrix->owns_blocks = true;
}


// End of file
