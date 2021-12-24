#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"

using namespace mfem;


void heat::LsqXtFemH1H1
:: setSpatialFeSpaceForHeatFlux(std::shared_ptr<Mesh> &spatialMesh)
{
    FiniteElementCollection *xH1Coll
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    auto spatialFeSpaceForHeatFlux
            = new FiniteElementSpace(spatialMesh.get(),
                                     xH1Coll, m_xDim);

    m_spatialFeCollections.Append(xH1Coll);
    m_spatialFeSpaces.Append(spatialFeSpaceForHeatFlux);
}

void heat::LsqXtFemH1H1
:: assembleSpatialMassForHeatFlux()
{
    BilinearForm *spatialMassForm
            = new BilinearForm(m_spatialFeSpaces[1]);
    spatialMassForm->AddDomainIntegrator
            (new VectorMassIntegrator);
    spatialMassForm->Assemble();
    spatialMassForm->Finalize();
    m_spatialMass2 = spatialMassForm->LoseMat();
    delete spatialMassForm;
}

void heat::LsqXtFemH1H1
:: assembleSpatialStiffnessForHeatFlux()
{
    BilinearForm *spatialStiffnessForm
            = new BilinearForm(m_spatialFeSpaces[1]);
    spatialStiffnessForm->AddDomainIntegrator
            (new heat::SpatialVectorStiffnessIntegrator);
    spatialStiffnessForm->Assemble();
    spatialStiffnessForm->Finalize();
    m_spatialStiffness2 = spatialStiffnessForm->LoseMat();
    delete spatialStiffnessForm;
}

void heat::LsqXtFemH1H1
:: assembleSpatialGradient()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    MixedBilinearForm *spatialGradientForm
            = new MixedBilinearForm(m_spatialFeSpaces[0],
                                    m_spatialFeSpaces[1]);
    spatialGradientForm->AddDomainIntegrator
            (new heat::SpatialVectorGradientIntegrator(&mediumCoeff));
    spatialGradientForm->Assemble();
    spatialGradientForm->Finalize();
    m_spatialGradient = spatialGradientForm->LoseMat();
    delete spatialGradientForm;
}

void heat::LsqXtFemH1H1
:: assembleSpatialDivergence()
{
    MixedBilinearForm *spatialDivergenceForm
            = new MixedBilinearForm(m_spatialFeSpaces[1],
                                    m_spatialFeSpaces[0]);
    spatialDivergenceForm->AddDomainIntegrator
            (new VectorDivergenceIntegrator);
    spatialDivergenceForm->Assemble();
    spatialDivergenceForm->Finalize();
    m_spatialDivergence = spatialDivergenceForm->LoseMat();
    delete spatialDivergenceForm;
}

void heat::LsqXtFemH1H1
:: assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis
(Vector& b, double t) const
{
    LinearForm source_form(m_spatialFeSpaces[1]);
    heat::SourceCoeff f(m_testCase);
    f.SetTime(t);
    source_form.AddDomainIntegrator
            (new heat::SpatialVectorDivergenceLFIntegrator(f,-1));
    source_form.Assemble();
    b = source_form.GetData();
}

// End of file
