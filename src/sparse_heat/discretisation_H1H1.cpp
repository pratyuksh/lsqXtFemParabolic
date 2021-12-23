#include "discretisation.hpp"

#include "assembly.hpp"
#include "../heat/coefficients.hpp"
#include "../heat/assembly.hpp"

using namespace mfem;


sparseHeat::LsqSparseXtFemH1H1
:: LsqSparseXtFemH1H1
(const nlohmann::json & config,
 std::shared_ptr<heat::TestCases> & testCase,
 const int numLevels,
 const int minTemporalLevel)
    : LsqSparseXtFem(config, testCase, numLevels, minTemporalLevel)
{}

sparseHeat::LsqSparseXtFemH1H1
:: LsqSparseXtFemH1H1
(const nlohmann::json & config,
 std::shared_ptr<heat::TestCases> & testCase,
 const int numLevels,
 const int minTemporalLevel,
 std::shared_ptr<mymfem::NestedMeshHierarchy> & spatialMeshHierachy)
    : LsqSparseXtFemH1H1(config, testCase, numLevels, minTemporalLevel)
{
    setNestedFEHierarchyAndSpatialBoundaryDofs(spatialMeshHierachy);
}

void sparseHeat::LsqSparseXtFemH1H1
:: setNestedFEHierarchyForHeatFlux
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    auto meshes = spatialNestedMeshHierarchy->getMeshes();
    assert((int) meshes.size() == m_numLevels);

    m_spatialNestedFEHierarchyHeatFlux
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
                (meshes[i].get(), xH1Coll, m_xDim);
        m_spatialNestedFEHierarchyHeatFlux->addFESpace(xH1Space);
    }
}

void sparseHeat::LsqSparseXtFemH1H1
:: assembleSpatialMassForHeatFlux()
{
    auto spatialVectorMassBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyHeatFlux);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialVectorMassIntegator
            = std::make_shared<sparseHeat::SpatialVectorMassIntegrator>();
    spatialVectorMassBilinearForm->addDomainIntegrator
            (spatialVectorMassIntegator);
    spatialVectorMassBilinearForm->assemble();
    m_spatialMass2 = spatialVectorMassBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1H1
:: assembleSpatialStiffnessForHeatFlux()
{
    auto spatialVectorStiffnessBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyHeatFlux);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialVectorStiffnessIntegator
            = std::make_shared<sparseHeat::SpatialVectorStiffnessIntegrator>();
    spatialVectorStiffnessBilinearForm->addDomainIntegrator
            (spatialVectorStiffnessIntegator);
    spatialVectorStiffnessBilinearForm->assemble();
    m_spatialStiffness2 = spatialVectorStiffnessBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1H1
:: assembleSpatialGradient()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    auto spatialVectorGradientBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (m_spatialNestedFEHierarchyTemperature,
             m_spatialNestedFEHierarchyHeatFlux);
    std::shared_ptr<mymfem::BlockMixedBilinearFormIntegrator>
            spatialVectorGradientIntegator
            = std::make_shared
            <sparseHeat::SpatialVectorGradientIntegrator>(mediumCoeff);
    spatialVectorGradientBilinearForm->addDomainIntegrator
            (spatialVectorGradientIntegator);
    spatialVectorGradientBilinearForm->assemble();
    m_spatialGradient = spatialVectorGradientBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1H1
:: assembleSpatialDivergence()
{
    auto spatialVectorDivergenceBilinearForm
            = std::make_unique<mymfem::BlockMixedBilinearForm>
            (m_spatialNestedFEHierarchyHeatFlux,
             m_spatialNestedFEHierarchyTemperature);
    std::shared_ptr<mymfem::BlockMixedBilinearFormIntegrator>
            spatialVectorDivergenceIntegator
            = std::make_shared<sparseHeat::SpatialVectorDivergenceIntegrator>();
    spatialVectorDivergenceBilinearForm->addDomainIntegrator
            (spatialVectorDivergenceIntegator);
    spatialVectorDivergenceBilinearForm->assemble();
    m_spatialDivergence = spatialVectorDivergenceBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFemH1H1
:: assembleSourceWithSpatialDivergenceOfHeatFluxBasisAtGivenTime
(double t,
 const std::shared_ptr<FiniteElementSpace>& spatialFes,
 Vector& b) const
{
    LinearForm sourceForm(spatialFes.get());
    heat::SourceCoeff sourceCoeff(m_testCase);
    sourceCoeff.SetTime(t);
    sourceForm.AddDomainIntegrator
            (new heat::SpatialVectorDivergenceLFIntegrator(sourceCoeff, -1));
    sourceForm.Assemble();
    b = sourceForm.GetData();
}

// End of file
