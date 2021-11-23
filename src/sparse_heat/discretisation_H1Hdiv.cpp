#include "discretisation.hpp"

#include "assembly.hpp"
#include "../heat/coefficients.hpp"


sparseHeat::LsqSparseXtFemH1Hdiv
:: LsqSparseXtFemH1Hdiv
(const nlohmann::json & config,
 std::shared_ptr<heat::TestCases> & testCase)
    : LsqSparseXtFem(config, testCase)
{}

sparseHeat::LsqSparseXtFemH1Hdiv
:: LsqSparseXtFemH1Hdiv
(const nlohmann::json & config,
 std::shared_ptr<heat::TestCases> & testCase,
 std::shared_ptr<mymfem::NestedMeshHierarchy> & spatialMeshHierachy)
    : LsqSparseXtFemH1Hdiv(config, testCase)
{
    setNestedFEHierarchyAndSpatialBoundaryDofs(spatialMeshHierachy);
}

void sparseHeat::LsqSparseXtFemH1Hdiv
:: setNestedFEHierarchyForHeatFlux
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    auto meshes = spatialNestedMeshHierarchy->getMeshes();
    assert((int) meshes.size() == m_numLevels);

    m_spatialNestedFEHierarchyHeatFlux
            = std::make_shared<mymfem::NestedFEHierarchy>
            (spatialNestedMeshHierarchy);

    m_xDim = meshes[0]->Dimension();
    FiniteElementCollection *xRTColl
            = new RT_FECollection(m_deg-1, m_xDim);
    m_xFecs.Append(xRTColl);

    for (int i=0; i<m_numLevels; i++)
    {
        auto xRTSpace = std::make_shared<FiniteElementSpace>
                (meshes[i].get(), xRTColl);
        m_spatialNestedFEHierarchyHeatFlux->addFESpace(xRTSpace);
    }
}

void sparseHeat::LsqSparseXtFemH1Hdiv
:: assembleSpatialMassForHeatFlux()
{
    auto spatialVectorFEMassBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyHeatFlux);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialVectorFEMassIntegator
            = std::make_shared<sparseHeat::SpatialVectorFEMassIntegrator>();
    spatialVectorFEMassBilinearForm->addDomainIntegrator
            (spatialVectorFEMassIntegator);
    spatialVectorFEMassBilinearForm->assemble();
    m_spatialMass2 = spatialVectorFEMassBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1Hdiv
:: assembleSpatialStiffnessForHeatFlux()
{
    auto spatialVectorFEStiffnessBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyHeatFlux);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialVectorFEStiffnessIntegator
            = std::make_shared<sparseHeat::SpatialVectorFEStiffnessIntegrator>();
    spatialVectorFEStiffnessBilinearForm->addDomainIntegrator
            (spatialVectorFEStiffnessIntegator);
    spatialVectorFEStiffnessBilinearForm->assemble();
    m_spatialStiffness2 = spatialVectorFEStiffnessBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1Hdiv
:: assembleSpatialGradient()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    auto spatialVectorFEGradientBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (m_spatialNestedFEHierarchyTemperature,
             m_spatialNestedFEHierarchyHeatFlux);
    std::shared_ptr<mymfem::BlockMixedBilinearFormIntegrator>
            spatialVectorFEGradientIntegator
            = std::make_shared
            <sparseHeat::SpatialVectorFEGradientIntegrator>(mediumCoeff);
    spatialVectorFEGradientBilinearForm->addDomainIntegrator
            (spatialVectorFEGradientIntegator);
    spatialVectorFEGradientBilinearForm->assemble();
    m_spatialGradient = spatialVectorFEGradientBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1Hdiv
:: assembleSpatialDivergence()
{
    auto spatialVectorFEDivergenceBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (m_spatialNestedFEHierarchyHeatFlux,
             m_spatialNestedFEHierarchyTemperature);
    std::shared_ptr<mymfem::BlockMixedBilinearFormIntegrator>
            spatialVectorFEDivergenceIntegator
            = std::make_shared<sparseHeat::SpatialVectorFEDivergenceIntegrator>();
    spatialVectorFEDivergenceBilinearForm->addDomainIntegrator
            (spatialVectorFEDivergenceIntegator);
    spatialVectorFEDivergenceBilinearForm->assemble();
    m_spatialDivergence = spatialVectorFEDivergenceBilinearForm->getBlockMatrix();
}

// End of file
