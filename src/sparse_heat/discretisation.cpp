#include "discretisation.hpp"

#include "assembly.hpp"
#include "utilities.hpp"
#include "../heat/coefficients.hpp"


sparseHeat::LsqSparseXtFem
:: LsqSparseXtFem (const nlohmann::json& config,
                   std::shared_ptr<heat::TestCases>& testCase)
    : m_config (config),
      m_testCase(testCase)
{
    m_deg = 1;
    if (config.contains("deg")) {
        m_deg = config["deg"];
        if (m_deg > 1) {
            std::cerr << "The sparse implementation has only been implemented "
                         "for polynomial degree 1" << std::endl;
            abort();
        }
    }

    m_endTime = 1.;
    if (config.contains("end_time")) {
        m_endTime = config["end_time"];
    }

    m_minTemporalLevel = 1;
    if (config.contains("min_temporal_level")) {
        m_minTemporalLevel = config["min_temporal_level"];
    }

    m_numLevels = 1;
    if (config.contains("num_levels")) {
        m_numLevels = config["num_levels"];
    }
    m_maxTemporalLevel = m_minTemporalLevel + m_numLevels - 1;
}

sparseHeat::LsqSparseXtFem
:: LsqSparseXtFem
(const nlohmann::json & config,
 std::shared_ptr<heat::TestCases> & testCase,
 std::shared_ptr<mymfem::NestedMeshHierarchy> & spatialMeshHierachy)
    : LsqSparseXtFem(config, testCase)
{
    setNestedFEHierarchyAndSpatialBoundaryDofs(spatialMeshHierachy);
}

sparseHeat::LsqSparseXtFem
:: ~LsqSparseXtFem()
{
    if (m_systemBlock11) {
        delete m_systemBlock11;
    }

    if (m_systemBlock12) {
        delete m_systemBlock12;
    }

    if (m_systemBlock21) {
        delete m_systemBlock21;
    }

    if (m_systemBlock22) {
        delete m_systemBlock22;
    }

    if (m_systemMatrix) {
        delete m_systemMatrix;
    }
}

void sparseHeat::LsqSparseXtFem
:: setNestedFEHierarchyAndSpatialBoundaryDofs
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    setNestedFEHierarchy(spatialNestedMeshHierarchy);
    setSpatialBoundaryDofs();
}

void sparseHeat::LsqSparseXtFem
:: setNestedFEHierarchy
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    setNestedFEHierarchyForTemperature(spatialNestedMeshHierarchy);
    setNestedFEHierarchyForHeatFlux(spatialNestedMeshHierarchy);
}

void sparseHeat::LsqSparseXtFem
:: setNestedFEHierarchyForTemperature
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    auto meshes = spatialNestedMeshHierarchy->getMeshes();
    assert((int) meshes.size() == m_numLevels);

    m_spatialNestedFEHierarchyTemperature
            = std::make_shared<mymfem::NestedFEHierarchy>
            (spatialNestedMeshHierarchy);

    m_xDim = meshes[0]->Dimension();
    FiniteElementCollection *xH1Coll
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    m_xFecs.Append(xH1Coll);

    for (int i=0; i<m_numLevels; i++)
    {
        auto xH1Space = std::make_shared<FiniteElementSpace>
                (meshes[i].get(), xH1Coll);
        m_spatialNestedFEHierarchyTemperature->addFESpace(xH1Space);
    }
}

void sparseHeat::LsqSparseXtFem
:: setSpatialBoundaryDofs()
{
    auto xMeshes
            = m_spatialNestedFEHierarchyTemperature->getMeshes();
    auto xFeSpaces
            = m_spatialNestedFEHierarchyTemperature->getFESpaces();

    auto tHierarchicalFESpacesSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto xFeSpacesSizes
            = m_spatialNestedFEHierarchyTemperature->getNumDims();
    auto sizeBuf = evalSpaceTimeBlockSizes(xFeSpacesSizes,
                                           tHierarchicalFESpacesSizes);
    auto offsetsBuf = evalBlockOffsets(sizeBuf);

    for (int m=0; m<m_numLevels; m++)
    {
        int i = getSpatialIndex(m);
        int shift = offsetsBuf[m];
        if (xMeshes[i]->bdr_attributes.Size())
        {
            Array<int> xEssBdrMarker(xMeshes[i]->bdr_attributes.Max());
            m_testCase->setBdryDirichlet(xEssBdrMarker);

            Array<int> xEssDofs;
            xFeSpaces[i]->GetEssentialTrueDofs(xEssBdrMarker, xEssDofs);

            int tNumDofs = tHierarchicalFESpacesSizes[m];
            int xNumDofs = xFeSpacesSizes[i];

            int xNumEssDofs = xEssDofs.Size();
            Array<int> myEssDofs(xNumEssDofs*tNumDofs);
            for (int j=0; j<tNumDofs; j++) {
                for (int k=0; k<xNumEssDofs; k++) {
                    myEssDofs[k + j*xNumEssDofs]
                            = xEssDofs[k] + j*xNumDofs + shift;
                }
            }
            m_essDofs.Append(myEssDofs);
//            myEssDofs.Print();
//            std::cout << "\n\n";
        }
    }
}

void sparseHeat::LsqSparseXtFem
:: assembleSystemSubMatrices()
{
    assembleTemporalInitial();
    assembleTemporalMass();
    assembleTemporalStiffness();
    assembleTemporalGradient();

    assembleSpatialMassForTemperature();
    assembleSpatialMassForHeatFlux();
    assembleSpatialStiffnessForTemperature();
    assembleSpatialStiffnessForHeatFlux();
    assembleSpatialGradient();
    assembleSpatialDivergence();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalInitial()
{
    auto temporalInitialAssembler
            = std::make_unique<sparseHeat::TemporalInitialMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalInitialAssembler->assemble();
    m_temporalInitial = temporalInitialAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalMass()
{
    auto temporalMassAssembler
            = std::make_unique<sparseHeat::TemporalMassMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalMassAssembler->assemble();
    m_temporalMass = temporalMassAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalStiffness()
{
    auto temporalStiffnessAssembler
            = std::make_unique<sparseHeat::TemporalStiffnessMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalStiffnessAssembler->assemble();
    m_temporalStiffness = temporalStiffnessAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalGradient()
{
    auto temporalGradientAssembler
            = std::make_unique<sparseHeat::TemporalGradientMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalGradientAssembler->assemble();
    m_temporalGradient = temporalGradientAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleSpatialMassForTemperature()
{
    auto spatialMassBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyTemperature);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialMassIntegator
            = std::make_shared<sparseHeat::SpatialMassIntegrator>();
    spatialMassBilinearForm->addDomainIntegrator(spatialMassIntegator);
    spatialMassBilinearForm->assemble();
    m_spatialMass1 = spatialMassBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleSpatialStiffnessForTemperature()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    auto spatialStiffnessBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyTemperature);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialStiffnessIntegator
            = std::make_shared<sparseHeat::SpatialStiffnessIntegrator>
            (mediumCoeff);
    spatialStiffnessBilinearForm->addDomainIntegrator
            (spatialStiffnessIntegator);
    spatialStiffnessBilinearForm->assemble();
    m_spatialStiffness1 = spatialStiffnessBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: buildSystemMatrix()
{
    buildSystemBlocks();

    Array<int> blockOffsets(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = m_systemBlock11->NumRows();
    blockOffsets[2] = m_systemBlock21->NumRows();
    blockOffsets.PartialSum();

    BlockMatrix buf(blockOffsets);
    buf.SetBlock(0, 0, m_systemBlock11);
    buf.SetBlock(0, 1, m_systemBlock12);
    buf.SetBlock(1, 0, m_systemBlock21);
    buf.SetBlock(1, 1, m_systemBlock22);

    m_systemMatrix = buf.CreateMonolithic();
    applyBCs(*m_systemMatrix);
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlocks()
{
    buildSystemBlock11();
    buildSystemBlock22();
    buildSystemBlock12();
    buildSystemBlock21();
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock11()
{
    auto temporalBlockSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto spatialBlockSizes
            = m_spatialNestedFEHierarchyTemperature->getNumDims();

    auto blockSizes = evalSpaceTimeBlockSizes(spatialBlockSizes,
                                              temporalBlockSizes);
    auto blockOffsets = evalBlockOffsets(blockSizes);

    auto buf = std::make_unique<BlockMatrix>(blockOffsets);
    // block11 = M^t_{\ell1, \ell2} \otimes K^{x;(1,1)}_{L-\ell1, L-\ell2}
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto mat = OuterProduct(m_temporalMass->GetBlock(i, j),
                                    m_spatialStiffness1->GetBlock(ii, jj));
            buf->SetBlock(i, j, mat);
        }
    }
    buf->owns_blocks = true;
    m_systemBlock11 = buf->CreateMonolithic();

    // block11 += (K^t_{\ell1, \ell2} + E^t_{\ell1, \ell2})
    //            \otimes M^{x; (1,1)}_{L-\ell1, L-\ell2}
    buf.reset(new BlockMatrix(blockOffsets));
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Add(m_temporalInitial->GetBlock(i, j),
                           m_temporalStiffness->GetBlock(i, j));
            auto mat = OuterProduct(*tmp, m_spatialMass1->GetBlock(ii, jj));
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    auto tmp = buf->CreateMonolithic();
    m_systemBlock11->Add(1., *tmp);
    delete tmp;
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock22()
{
    auto temporalBlockSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto spatialBlockSizes
            = m_spatialNestedFEHierarchyHeatFlux->getNumDims();

    auto blockSizes = evalSpaceTimeBlockSizes(spatialBlockSizes,
                                              temporalBlockSizes);
    auto blockOffsets = evalBlockOffsets(blockSizes);

    auto buf = std::make_unique<BlockMatrix>(blockOffsets);
    // block22 = M^t_{\ell1, \ell2}
    // \otimes (M^{x; (2,2)}_{L-\ell1, L-\ell2}
    //        + K^{x; (2,2)}_{L-\ell1, L-\ell2})
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Add(m_spatialMass2->GetBlock(ii, jj),
                           m_spatialStiffness2->GetBlock(ii, jj));
            auto mat = OuterProduct(m_temporalMass->GetBlock(i, j), *tmp);
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    m_systemBlock22 = buf->CreateMonolithic();
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock12()
{
    auto temporalBlockSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto spatialRowBlockSizes
            = m_spatialNestedFEHierarchyTemperature->getNumDims();
    auto spatialColBlockSizes
            = m_spatialNestedFEHierarchyHeatFlux->getNumDims();

    auto rowBlockSizes = evalSpaceTimeBlockSizes(spatialRowBlockSizes,
                                                 temporalBlockSizes);
    auto colBlockSizes = evalSpaceTimeBlockSizes(spatialColBlockSizes,
                                                 temporalBlockSizes);

    auto rowBlockOffsets = evalBlockOffsets(rowBlockSizes);
    auto colBlockOffsets = evalBlockOffsets(colBlockSizes);

//    rowBlockSizes.Print();
//    colBlockSizes.Print();

    auto buf = std::make_unique<BlockMatrix>(rowBlockOffsets,
                                             colBlockOffsets);
    // block12 = (C^t_{\ell1, \ell2})^{\top}
    //           \otimes B^{x; (1,2)}_{L-\ell1, L-\ell2})
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Transpose(m_temporalGradient->GetBlock(j, i));
            auto mat = OuterProduct
                    (*tmp, m_spatialDivergence->GetBlock(ii, jj));
//            std::cout << i << "\t"
//                      << j << "\t"
//                      << mat->NumRows() << "\t"
//                      << mat->NumCols() << "\t"
//                      << rowBlockSizes[i] << "\t"
//                      << colBlockSizes[j] << std::endl;
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    m_systemBlock12 = buf->CreateMonolithic();

    // block12 += (M^t_{\ell1, \ell2})^{\top}
    //           \otimes C^{x; (1,2)}_{L-\ell1, L-\ell2}
    buf.reset(new BlockMatrix(rowBlockOffsets, colBlockOffsets));
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Transpose(m_spatialGradient->GetBlock(jj, ii));
            auto mat = OuterProduct(m_temporalMass->GetBlock(i, j), *tmp);
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    auto tmp = buf->CreateMonolithic();
    m_systemBlock12->Add(1., *tmp);
    delete tmp;
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock21()
{
    m_systemBlock21 = Transpose(*m_systemBlock12);
}

int sparseHeat::LsqSparseXtFem
:: getSpatialIndex(int i) const {
    return m_numLevels - i - 1;
}

void sparseHeat::LsqSparseXtFem
:: assembleRhs(std::shared_ptr<BlockVector>& B) const
{
    (*B) = 0.;

    Vector &B0 = B->GetBlock(0);
    assembleICs(B0);

//    Vector &B1 = B->GetBlock(1);
//    assembleSource(B1);

    applyBCs(*B);
}

void sparseHeat::LsqSparseXtFem
:: assembleICs(Vector& b) const
{
    auto xFeSpaces
            = m_spatialNestedFEHierarchyTemperature->getFESpaces();
    int coarsestTemporalIndex = 0;
    int finestSpatialIndex = getSpatialIndex(coarsestTemporalIndex);

    LinearForm icsForm(xFeSpaces[finestSpatialIndex].get());
    heat::InitialTemperatureCoeff u0(m_testCase);
    icsForm.AddDomainIntegrator
            (new DomainLFIntegrator(u0));
    icsForm.Assemble();
    Vector buf(icsForm.GetData(),
               xFeSpaces[finestSpatialIndex]->GetTrueVSize());

    b.SetVector(buf, 0);
}

void sparseHeat::LsqSparseXtFem
:: assembleSource(Vector& b) const
{}

void sparseHeat::LsqSparseXtFem
:: applyBCs(SparseMatrix &A) const
{
    for (int k=0; k<m_essDofs.Size(); k++) {
        A.EliminateRowCol(m_essDofs[k]);
    }
}

void sparseHeat::LsqSparseXtFem
:: applyBCs(BlockVector &B) const
{
    auto& b = B.GetBlock(0);
    b.SetSubVector(m_essDofs, 0.);
}

// End of file
