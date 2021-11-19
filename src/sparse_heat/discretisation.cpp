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

// Destructor
sparseHeat::LsqSparseXtFem
:: ~LsqSparseXtFem()
{
    // reset();
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
    // block11 = M_t \otimes K_x^{1,1}
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
    m_block11 = buf->CreateMonolithic();

    // block11 += (K_t + E_t) \otimes M_x^{1,1}
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
    m_block11->Add(1., *tmp);
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
    // block22 = M_t \otimes (M_x^{2,2} + K_x^{2,2})
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
    m_block22 = buf->CreateMonolithic();
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

    rowBlockSizes.Print();
    colBlockSizes.Print();

    auto buf = std::make_unique<BlockMatrix>(rowBlockOffsets,
                                             colBlockOffsets);
    // block12 = C_t^{\top} \otimes B_x^{1,2}
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
    m_block12 = buf->CreateMonolithic();

    // block12 += M_t \otimes C_x^{1,2}
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
    m_block12->Add(1., *tmp);
    delete tmp;
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock21()
{
    m_block21 = Transpose(*m_block12);
}

// End of file
