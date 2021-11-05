#include "my_bilinearForms.hpp"
#include "utilities.hpp"


// //! Assembles the domain integrators, interior face integrators
// //! and boundary face integrators
//void mymfem::MyBilinearForm :: MyAssemble (int skip_zeros)
//{
//    // allocate memory
//    if (mat == nullptr) {
//        AllocMat();
//    }

//    // assemble domain integrators
//    if (mydbfi.Size())
//    {
//        DenseMatrix elmat;
//        ElementTransformation *eltrans = nullptr;

//        for (int i = 0; i < fes -> GetNE(); i++)
//        {
//            fes->GetElementVDofs(i, vdofs);

//            const FiniteElement &fe = *fes->GetFE(i);
//            eltrans = fes->GetElementTransformation(i);
//            for (int k = 0; k < mydbfi.Size(); k++)
//            {
//                mydbfi[k]->AssembleElementMatrix(fe, *eltrans, elmat);
//                mat->AddSubMatrix(vdofs, vdofs, elmat, skip_zeros);
//            }
//        }
//    }

//    // assemble interior face integrators
//    if (myifbfi.Size())
//    {
//        DenseMatrix elmat;
//        FaceElementTransformations *ftr = nullptr;
//        Array<int> vdofs2;
//        Mesh *mesh = fes -> GetMesh();

//        int nfaces = mesh->GetNumFaces();
//        for (int i = 0; i < nfaces; i++)
//        {
//            ftr = mesh -> GetInteriorFaceTransformations (i);
//            if (ftr != nullptr)
//            {
//                const FiniteElement &fe1 = *(fes->GetFE (ftr->Elem1No));
//                const FiniteElement &fe2 = *(fes->GetFE (ftr->Elem2No));

//                fes -> GetElementVDofs (ftr -> Elem1No, vdofs);
//                fes -> GetElementVDofs (ftr -> Elem2No, vdofs2);
//                vdofs.Append (vdofs2);
//                for (int k = 0; k < myifbfi.Size(); k++)
//                {
//                    myifbfi[k] -> AssembleFaceMatrix
//                            (fe1, fe2, *ftr, elmat);

//                    mat -> AddSubMatrix (vdofs, vdofs, elmat, skip_zeros);
//                }
//            }
//        }
//    }

//    // assemble boundary face integrators
//    if (mybfbfi.Size())
//    {
//        DenseMatrix elmat;
//        FaceElementTransformations *ftr = nullptr;
//        Mesh *mesh = fes -> GetMesh();

//        int nbfaces = mesh->GetNBE();
//        for (int i = 0; i < nbfaces; i++)
//        {
//            ftr = mesh -> GetBdrFaceTransformations(i);
//            const FiniteElement &fe = *(fes->GetFE (ftr->Elem1No));

//            fes -> GetElementVDofs (ftr -> Elem1No, vdofs);
//            for (int k = 0; k < mybfbfi.Size(); k++)
//            {
//                mybfbfi[k] -> AssembleFaceMatrix
//                        (fe, fe, *ftr, elmat);
//                mat -> AddSubMatrix (vdofs, vdofs, elmat, skip_zeros);
//            }
//        }
//    }
//}


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
:: assemble (int skip_zeros)
{
    // if not initialized, allocate block matrix
    if (!m_blockMatrix) {
        this->allocateBlockMatrix();
    }

    if (mydbfi.size())
    {
        auto feSpaces = m_nestedFEHierarchy->getFESpaces();
        auto meshes = m_nestedFEHierarchy->getMeshes();
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
                for (int k = 0; k < mydbfi.size(); k++)
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

                auto& mat = m_blockMatrix->GetBlock(m,n);
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
                        for (int k = 0; k < mydbfi.size(); k++)
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
        m_blockMatrix->Finalize(skip_zeros, true);

        // upper-triangle blocks
        for (int n=1; n<m_numLevels; n++)
            for (int m=0; m<n; m++)
            {
                auto matT = Transpose(m_blockMatrix->GetBlock(n,m));
                m_blockMatrix->GetBlock(m,n)
                        = *matT;
                delete matT;
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


// End of file
