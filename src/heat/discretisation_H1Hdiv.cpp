#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"

using namespace mfem;


void heat::LsqXtFemH1Hdiv
:: setSpatialFeSpaceForHeatFlux(std::shared_ptr<Mesh> &spatialMesh)
{
    FiniteElementCollection *xRTColl
            = new RT_FECollection(m_deg-1, m_xDim);
    auto spatialFeSpaceForHeatFlux
            = new FiniteElementSpace(spatialMesh.get(),
                                     xRTColl);

    m_spatialFeCollections.Append(xRTColl);
    m_spatialFeSpaces.Append(spatialFeSpaceForHeatFlux);
}

void heat::LsqXtFemH1Hdiv
:: assembleSpatialMassForHeatFlux()
{
    BilinearForm *spatialMassForm
            = new BilinearForm(m_spatialFeSpaces[1]);
    spatialMassForm->AddDomainIntegrator
            (new VectorFEMassIntegrator);
    spatialMassForm->Assemble();
    spatialMassForm->Finalize();
    m_spatialMass2 = spatialMassForm->LoseMat();
    delete spatialMassForm;
}

void heat::LsqXtFemH1Hdiv
:: assembleSpatialStiffnessForHeatFlux()
{
    BilinearForm *spatialStiffnessForm
            = new BilinearForm(m_spatialFeSpaces[1]);
    spatialStiffnessForm->AddDomainIntegrator
            (new heat::SpatialVectorFEStiffnessIntegrator);
    spatialStiffnessForm->Assemble();
    spatialStiffnessForm->Finalize();
    m_spatialStiffness2 = spatialStiffnessForm->LoseMat();
    delete spatialStiffnessForm;
}

void heat::LsqXtFemH1Hdiv
:: assembleSpatialGradient()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    MixedBilinearForm *spatialGradientForm
            = new MixedBilinearForm(m_spatialFeSpaces[0],
                                    m_spatialFeSpaces[1]);
    spatialGradientForm->AddDomainIntegrator
            (new heat::SpatialVectorFEGradientIntegrator(&mediumCoeff));
    spatialGradientForm->Assemble();
    spatialGradientForm->Finalize();
    m_spatialGradient = spatialGradientForm->LoseMat();
    delete spatialGradientForm;
}

void heat::LsqXtFemH1Hdiv
:: assembleSpatialDivergence()
{
    MixedBilinearForm *spatialDivergenceForm
            = new MixedBilinearForm(m_spatialFeSpaces[1],
                                    m_spatialFeSpaces[0]);
    spatialDivergenceForm->AddDomainIntegrator
            (new VectorFEDivergenceIntegrator);
    spatialDivergenceForm->Assemble();
    spatialDivergenceForm->Finalize();
    m_spatialDivergence = spatialDivergenceForm->LoseMat();
    delete spatialDivergenceForm;
}

void heat::LsqXtFemH1Hdiv
:: assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis
(Vector& b, double t) const
{
    LinearForm source_form(m_spatialFeSpaces[1]);
    heat::SourceCoeff f(m_testCase);
    f.SetTime(t);
    source_form.AddDomainIntegrator
            (new heat::SpatialVectorFEDivergenceLFIntegrator(f,-1));
    source_form.Assemble();
    b = source_form.GetData();
}

// End of file
