#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"

using namespace mfem;


// Sets the FE space and boundary conditions
/*void heat::LsqXtFemH1Hdiv
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
    FiniteElementCollection *xH1_coll
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    FiniteElementSpace *xV_space
            = new FiniteElementSpace(xMesh.get(), xH1_coll);
    
    FiniteElementCollection *xHdiv_coll
            = new RT_FECollection(m_deg-1, m_xDim);
    FiniteElementSpace *xR_space
            = new FiniteElementSpace(xMesh.get(),
                                     xHdiv_coll);
    
    int tdimV = m_temporalFeSpace->GetTrueVSize();
    int xdimV = xV_space->GetTrueVSize();
    int xdimR = xR_space->GetTrueVSize();
#ifndef NDEBUG
    std::cout << "\n\tNumber of degrees of freedom "
                 "in tV_space: "
         << tdimV << std::endl;
    std::cout << "\tNumber of degrees of freedom "
                 "in xV_space: "
         << xdimV << std::endl;
    std::cout << "\tNumber of degrees of freedom"
                 " in xR_space: "
         << xdimR << std::endl;
#endif
    m_spatialFeCollections.Append(xH1_coll);
    m_spatialFeCollections.Append(xHdiv_coll);
    
    m_spatialFeSpaces.Append(xV_space);
    m_spatialFeSpaces.Append(xR_space);

    // Define the two BlockStructure of the problem
    m_blockOffsets.SetSize(3);
    m_blockOffsets[0] = 0;
    m_blockOffsets[1] = xdimV*tdimV;
    m_blockOffsets[2] = xdimR*tdimV;
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

        xV_space->GetEssentialTrueDofs(m_spatialEssentialBoundaryMarker,
                                       m_spatialEssentialDofs);

        int xNEssDofs = m_spatialEssentialDofs.Size();
        m_essentialDofs.SetSize(xNEssDofs*tdimV);
        for (int j=0; j<tdimV; j++)
            for (int i=0; i<xNEssDofs; i++)
                m_essentialDofs[i + j*xNEssDofs]
                        = j*xdimV + m_spatialEssentialDofs[i];
    }
}*/

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

// Assembles medium-independent system matrices
/*void heat::LsqXtFemH1Hdiv
:: assembleSystemMediumIndependent()
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

    // xMass1 form, xV_space
    BilinearForm *xMass1_form
            = new BilinearForm(m_spatialFeSpaces[0]);
    xMass1_form->AddDomainIntegrator(new MassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_spatialMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xMass2 form, xR_space
    BilinearForm *xMass2_form
            = new BilinearForm(m_spatialFeSpaces[1]);
    xMass2_form->AddDomainIntegrator
            (new VectorFEMassIntegrator);
    xMass2_form->Assemble();
    xMass2_form->Finalize();
    m_spatialMass2 = xMass2_form->LoseMat();
    delete xMass2_form;

    // xStiff2 form, xR_space
    BilinearForm *xStiff2_form
            = new BilinearForm(m_spatialFeSpaces[1]);
    xStiff2_form->AddDomainIntegrator
            (new heat::SpatialVectorFEStiffnessIntegrator);
    xStiff2_form->Assemble();
    xStiff2_form->Finalize();
    m_spatialStiffness2 = xStiff2_form->LoseMat();
    delete xStiff2_form;

    // xDiv form
    MixedBilinearForm *xDiv_form
            = new MixedBilinearForm(m_spatialFeSpaces[1],
                                    m_spatialFeSpaces[0]);
    xDiv_form->AddDomainIntegrator
            (new VectorFEDivergenceIntegrator);
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
void heat::LsqXtFemH1Hdiv :: assembleSystemMediumDependent()
{
    // assemble matrices in space

    heat::MediumTensorCoeff medCoeff(m_testCase);

    // xStiff1 form, xV_space
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
            (new heat::SpatialVectorFEGradientIntegrator(&medCoeff));
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

void heat::LsqXtFemH1Hdiv
:: assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis(Vector& b, double t) const
{
    LinearForm source_form(m_spatialFeSpaces[1]);
    heat::SourceCoeff f(m_testCase);
    f.SetTime(t);
    source_form.AddDomainIntegrator
            (new heat::SpatialVectorFEDivergenceLFIntegrator(f,-1));
    source_form.Assemble();
    b = source_form.GetData();
}

// Assembles matrices for L2-projection
/*void heat::LsqXtFemH1Hdiv :: assembleProjector()
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

    // xMass1 form, xV_space
    BilinearForm *xMass1_form
            = new BilinearForm(m_spatialFeSpaces[0]);
    xMass1_form->AddDomainIntegrator(new MassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_spatialMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xMass2 form, xR_space
    BilinearForm *xMass2_form
            = new BilinearForm(m_spatialFeSpaces[1]);
    xMass2_form->AddDomainIntegrator
            (new VectorFEMassIntegrator);
    xMass2_form->Assemble();
    xMass2_form->Finalize();
    m_spatialMass2 = xMass2_form->LoseMat();
    delete xMass2_form;

    //compute blocks
    m_systemBlock11 = OuterProduct(*m_temporalMass, *m_spatialMass1);
    m_systemBlock11 = OuterProduct(*m_temporalMass, *m_spatialMass2);
}

void heat::LsqXtFemH1Hdiv
:: assembleSpatialProjectionOfRhs(Vector& b1, Vector& b2,
                            double t) const
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
            (new VectorFEDomainLFIntegrator(qCoeff));
    qSol_form.Assemble();
    b2 = qSol_form.GetData();
}*/

// End of file
