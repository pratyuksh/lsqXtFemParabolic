#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"

using namespace mfem;


// Sets the FE space and boundary conditions
/*void heat::LsqXtFemH1H1
:: setFESpacesAndSpatialBoundaryDofs(std::shared_ptr<Mesh>& tMesh,
                                     std::shared_ptr<Mesh>& xMesh)
{
    // reset initialized variables
    reset();

    m_temporalFeCollection = new H1_FECollection(m_deg, 1,
                                 BasisType::GaussLobatto);
    m_temporalFeSpace
            = new FiniteElementSpace(tMesh.get(), m_temporalFeCollection);
    
    m_xDim = xMesh->Dimension();
    FiniteElementCollection *xH1_coll1
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    FiniteElementSpace *xV1space
            = new FiniteElementSpace(xMesh.get(), xH1_coll1);
    
    FiniteElementCollection *xH1_coll2
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    FiniteElementSpace *xV2space
            = new FiniteElementSpace(xMesh.get(),
                                     xH1_coll2, m_xDim);
    
    int tdimV = m_temporalFeSpace->GetTrueVSize();
    int xdimV1 = xV1space->GetTrueVSize();
    int xdimV2 = xV2space->GetTrueVSize();
#ifndef NDEBUG
    std::cout << "\n\tNumber of degrees of freedom "
                 "in tV_space: "
         << tdimV << std::endl;
    std::cout << "\tNumber of degrees of freedom "
                 "in xV1space: "
         << xdimV1 << std::endl;
    std::cout << "\tNumber of degrees of freedom"
                 " in xV2space: "
         << xdimV2 << std::endl;
#endif
    m_spatialFeCollections.Append(xH1_coll1);
    m_spatialFeCollections.Append(xH1_coll2);
    
    m_spatialFeSpaces.Append(xV1space);
    m_spatialFeSpaces.Append(xV2space);

    // Define the two BlockStructure of the problem.
    m_blockOffsets.SetSize(3);
    m_blockOffsets[0] = 0;
    m_blockOffsets[1] = xdimV1*tdimV;
    m_blockOffsets[2] = xdimV2*tdimV;
    m_blockOffsets.PartialSum();
#ifndef NDEBUG
    std::cout << "\tBlock offsets:" << "\t"
              << m_blockOffsets[0] << "\t"
              << m_blockOffsets[1] << "\t"
              << m_blockOffsets[2] << std::endl;
#endif
    // Mark boundary dofs for xV_space
    if (xMesh->bdr_attributes.Size())
    {
        m_spatialEssentialBoundaryMarker.SetSize(xMesh->bdr_attributes.Max());
        m_testCase->setBdryDirichlet(m_spatialEssentialBoundaryMarker);

        xV1space->GetEssentialTrueDofs(m_spatialEssentialBoundaryMarker,
                                       m_spatialEssentialDofs);

        int xNEssDofs = m_spatialEssentialDofs.Size();
        m_essentialDofs.SetSize(xNEssDofs*tdimV);
        for (int j=0; j<tdimV; j++)
            for (int i=0; i<xNEssDofs; i++)
                m_essentialDofs[i + j*xNEssDofs]
                        = j*xdimV1 + m_spatialEssentialDofs[i];
    }
}*/

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

// Assembles medium-independent system matrices
/*void heat::LsqXtFemH1H1 :: assembleSystemMediumIndependent()
{
    // assemble matrices in time

    // tMass form
    BilinearForm *tMass_form = new BilinearForm(m_temporalFeSpace);
    tMass_form->AddDomainIntegrator(new MassIntegrator);
    tMass_form->Assemble();
    tMass_form->Finalize();
    m_temporalMass = tMass_form->LoseMat();
    delete tMass_form;

    // tStiffness form
    BilinearForm *tStiff_form = new BilinearForm(m_temporalFeSpace);
    tStiff_form->AddDomainIntegrator(new DiffusionIntegrator);
    tStiff_form->Assemble();
    tStiff_form->Finalize();
    m_temporalStiffness = tStiff_form->LoseMat();
    delete tStiff_form;

    // tGrad form
    BilinearForm *tGrad_form = new BilinearForm(m_temporalFeSpace);
    tGrad_form->AddDomainIntegrator
            (new heat::TemporalGradientIntegrator);
    tGrad_form->Assemble();
    tGrad_form->Finalize();
    m_temporalGradient = tGrad_form->LoseMat();
    delete tGrad_form;

    // assemble matrices in space

    // xMass1 form, xV1space
    BilinearForm *xMass1_form
            = new BilinearForm(m_spatialFeSpaces[0]);
    xMass1_form->AddDomainIntegrator(new MassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_spatialMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xMass2 form, xV2space
    BilinearForm *xMass2_form
            = new BilinearForm(m_spatialFeSpaces[1]);
    xMass2_form->AddDomainIntegrator
            (new VectorMassIntegrator);
    xMass2_form->Assemble();
    xMass2_form->Finalize();
    m_spatialMass2 = xMass2_form->LoseMat();
    delete xMass2_form;

    // xStiff2 form, xV2space
    BilinearForm *xStiff2_form
            = new BilinearForm(m_spatialFeSpaces[1]);
    xStiff2_form->AddDomainIntegrator
            (new heat::SpatialVectorStiffnessIntegrator);
    xStiff2_form->Assemble();
    xStiff2_form->Finalize();
    m_spatialStiffness2 = xStiff2_form->LoseMat();
    delete xStiff2_form;

    // xDiv form
    MixedBilinearForm *xDiv_form
            = new MixedBilinearForm(m_spatialFeSpaces[1],
                                    m_spatialFeSpaces[0]);
    xDiv_form->AddDomainIntegrator
            (new VectorDivergenceIntegrator);
    xDiv_form->Assemble();
    xDiv_form->Finalize();
    m_spatialDivergence = xDiv_form->LoseMat();
    delete xDiv_form;

    // matrix for initial conditions
    Vector tCanonicalBasis(m_temporalFeSpace->GetVSize());
    tCanonicalBasis = 0.0;
    tCanonicalBasis(0) = 1;
    SparseMatrix *tInit = new SparseMatrix(tCanonicalBasis);

    // compute blocks
    SparseMatrix *tmpBuf;

    // block at row 1, column 1
    tmpBuf = Add(*m_spatialMass2, *m_spatialStiffness2);
    m_systemBlock22 = OuterProduct(*m_temporalMass, *tmpBuf);
    delete tmpBuf;
    clear(m_spatialMass2);
    clear(m_spatialStiffness2);

    // block at row 0, column 0
    tmpBuf = Add(*tInit, *m_temporalStiffness);
    m_materialIndependentSystemBlock11 = OuterProduct(*tmpBuf, *m_spatialMass1);
    delete tmpBuf;
    clear(tInit);
    clear(m_temporalStiffness);
    clear(m_spatialMass1);

    // block at row 0, column 1
    tmpBuf = Transpose(*m_temporalGradient);
    m_materialIndependentSystemBlock12 = OuterProduct(*tmpBuf, *m_spatialDivergence);
    delete tmpBuf;
    clear(m_temporalGradient);
    clear(m_spatialDivergence);
}

// Assembles medium-dependent system matrices
void heat::LsqXtFemH1H1 :: assembleSystemMediumDependent()
{
    // assemble matrices in space

    heat::MediumTensorCoeff medCoeff(m_testCase);

    // xStiff1 form, xV1space
    BilinearForm *xStiff1_form
            = new BilinearForm(m_spatialFeSpaces[0]);
    xStiff1_form->AddDomainIntegrator
            (new heat::SpatialStiffnessIntegrator(&medCoeff));
    xStiff1_form->Assemble();
    xStiff1_form->Finalize();
    m_spatialStiffness1 = xStiff1_form->LoseMat();
    delete xStiff1_form;

    // xGrad form
    MixedBilinearForm *xGrad_form
            = new MixedBilinearForm(m_spatialFeSpaces[0],
                                    m_spatialFeSpaces[1]);
    xGrad_form->AddDomainIntegrator
            (new heat::SpatialVectorGradientIntegrator(&medCoeff));
    xGrad_form->Assemble();
    xGrad_form->Finalize();
    m_spatialGradient = xGrad_form->LoseMat();
    delete xGrad_form;

    // compute blocks
    SparseMatrix *tmp1, *tmp2, *tmpBuf;

    // block at row 0, column 0
    clear(m_systemBlock11);
    tmpBuf = OuterProduct(*m_temporalMass, *m_spatialStiffness1);
    m_systemBlock11 = Add(*m_materialIndependentSystemBlock11, *tmpBuf);
    delete tmpBuf;
    clear(m_spatialStiffness1);

    // block at row 0, column 1
    clear(m_systemBlock12);
    tmp2 = Transpose(*m_spatialGradient);
    tmp1 = OuterProduct(*m_temporalMass, *tmp2);
    m_systemBlock12 = Add(-1, *tmp1, -1, *m_materialIndependentSystemBlock12);
    delete tmp1;
    delete tmp2;
    clear(m_spatialGradient);

    // block at row 1, column 0
    clear(m_systemBlock21);
    m_systemBlock21 = Transpose(*m_systemBlock12);
}*/

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

// Assembles matrices for L2-projection
/*void heat::LsqXtFemH1H1 :: assembleProjector()
{
    // assemble matrices in time

    // tMass form
    BilinearForm *tMass_form = new BilinearForm(m_temporalFeSpace);
    tMass_form->AddDomainIntegrator(new MassIntegrator);
    tMass_form->Assemble();
    tMass_form->Finalize();
    m_temporalMass = tMass_form->LoseMat();
    delete tMass_form;

    // assemble matrices in space

    // xMass1 form, xV1space
    BilinearForm *xMass1_form
            = new BilinearForm(m_spatialFeSpaces[0]);
    xMass1_form->AddDomainIntegrator(new MassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_spatialMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xMass2 form, xV2space
    BilinearForm *xMass2_form
            = new BilinearForm(m_spatialFeSpaces[1]);
    xMass2_form->AddDomainIntegrator
            (new VectorMassIntegrator);
    xMass2_form->Assemble();
    xMass2_form->Finalize();
    m_spatialMass2 = xMass2_form->LoseMat();
    delete xMass2_form;

    // compute blocks
    m_systemBlock11 = OuterProduct(*m_temporalMass, *m_spatialMass1);
    m_systemBlock11 = OuterProduct(*m_temporalMass, *m_spatialMass2);
}

void heat::LsqXtFemH1H1
:: assembleSpatialProjectionOfRhs(Vector& b1, Vector& b2, double t) const
{
    LinearForm uSol_form(m_spatialFeSpaces[0]);
    heat::ExactTemperatureCoeff uCoeff(m_testCase);
    uCoeff.SetTime(t);
    uSol_form.AddDomainIntegrator
            (new DomainLFIntegrator(uCoeff));
    uSol_form.Assemble();
    b1 = uSol_form.GetData();

    LinearForm qSol_form(m_spatialFeSpaces[1]);
    heat::ExactHeatFluxCoeff qCoeff(m_testCase);
    qCoeff.SetTime(t);
    qSol_form.AddDomainIntegrator
            (new VectorDomainLFIntegrator(qCoeff));
    qSol_form.Assemble();
    b2 = qSol_form.GetData();
}*/

// End of file
