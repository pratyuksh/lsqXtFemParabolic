#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"


// Sets the FE space and boundary conditions
void heat::LsqXtFemH1Hdiv :: set(std::shared_ptr<Mesh>& tMesh,
                                 std::shared_ptr<Mesh>& xMesh)
{
    // reset initialized variables
    reset();

    m_tFec = new H1_FECollection(m_deg, 1,
                                 BasisType::GaussLobatto);
    m_tFespace
            = new FiniteElementSpace(tMesh.get(), m_tFec);
    
    m_xndim = xMesh->Dimension();
    FiniteElementCollection *xH1_coll
            = new H1_FECollection(m_deg, m_xndim,
                                  BasisType::GaussLobatto);
    FiniteElementSpace *xV_space
            = new FiniteElementSpace(xMesh.get(), xH1_coll);
    
    FiniteElementCollection *xHdiv_coll
            = new RT_FECollection(m_deg-1, m_xndim);
    FiniteElementSpace *xR_space
            = new FiniteElementSpace(xMesh.get(),
                                     xHdiv_coll);
    
    int tdimV = m_tFespace->GetTrueVSize();
    int xdimV = xV_space->GetTrueVSize();
    int xdimR = xR_space->GetTrueVSize();
#ifdef MYVERBOSE
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
    m_xFecs.Append(xH1_coll);
    m_xFecs.Append(xHdiv_coll);
    
    m_xFespaces.Append(xV_space);
    m_xFespaces.Append(xR_space);

    // Define the two BlockStructure of the problem
    m_block_offsets.SetSize(3);
    m_block_offsets[0] = 0;
    m_block_offsets[1] = xdimV*tdimV;
    m_block_offsets[2] = xdimR*tdimV;
    m_block_offsets.PartialSum();
#ifdef MYVERBOSE
    std::cout << "\tBlock offsets:" << "\t"
              << m_block_offsets[0] << "\t"
              << m_block_offsets[1] << "\t"
              << m_block_offsets[2] << std::endl;
#endif
    // Mark boundary dofs for xV_space
    if (xMesh->bdr_attributes.Size())
    {
        m_xEssBdrMarker.SetSize(xMesh->bdr_attributes.Max());
        m_testCase->setBdryDirichlet(m_xEssBdrMarker);

        xV_space->GetEssentialTrueDofs(m_xEssBdrMarker,
                                       m_xEssTdofList);

        int xNEssDofs = m_xEssTdofList.Size();
        m_essTdofList.SetSize(xNEssDofs*tdimV);
        for (int j=0; j<tdimV; j++)
            for (int i=0; i<xNEssDofs; i++)
                m_essTdofList[i + j*xNEssDofs]
                        = j*xdimV + m_xEssTdofList[i];
    }
}

// Assembles medium-independent system matrices
void heat::LsqXtFemH1Hdiv
:: assembleSystemMediumIndependent()
{
    // assemble matrices in time

    // tMass form
    BilinearForm *tMass_form = new BilinearForm(m_tFespace);
    tMass_form->AddDomainIntegrator(new MassIntegrator);
    tMass_form->Assemble();
    tMass_form->Finalize();
    m_tMass = tMass_form->LoseMat();
    delete tMass_form;

    // tStiffness form
    BilinearForm *tStiff_form = new BilinearForm(m_tFespace);
    tStiff_form->AddDomainIntegrator(new DiffusionIntegrator);
    tStiff_form->Assemble();
    tStiff_form->Finalize();
    m_tStiff = tStiff_form->LoseMat();
    delete tStiff_form;

    // tGrad form
    BilinearForm *tGrad_form = new BilinearForm(m_tFespace);
    tGrad_form->AddDomainIntegrator
            (new heat::GradientIntegrator);
    tGrad_form->Assemble();
    tGrad_form->Finalize();
    m_tGrad = tGrad_form->LoseMat();
    delete tGrad_form;

    // assemble matrices in space

    // xMass1 form, xV_space
    BilinearForm *xMass1_form
            = new BilinearForm(m_xFespaces[0]);
    xMass1_form->AddDomainIntegrator(new MassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_xMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xMass2 form, xR_space
    BilinearForm *xMass2_form
            = new BilinearForm(m_xFespaces[1]);
    xMass2_form->AddDomainIntegrator
            (new VectorFEMassIntegrator);
    xMass2_form->Assemble();
    xMass2_form->Finalize();
    m_xMass2 = xMass2_form->LoseMat();
    delete xMass2_form;

    // xStiff2 form, xR_space
    BilinearForm *xStiff2_form
            = new BilinearForm(m_xFespaces[1]);
    xStiff2_form->AddDomainIntegrator
            (new heat::VectorFEStiffnessIntegrator);
    xStiff2_form->Assemble();
    xStiff2_form->Finalize();
    m_xStiff2 = xStiff2_form->LoseMat();
    delete xStiff2_form;

    // xDiv form
    MixedBilinearForm *xDiv_form
            = new MixedBilinearForm(m_xFespaces[1],
                                    m_xFespaces[0]);
    xDiv_form->AddDomainIntegrator
            (new VectorFEDivergenceIntegrator);
    xDiv_form->Assemble();
    xDiv_form->Finalize();
    m_xDiv = xDiv_form->LoseMat();
    delete xDiv_form;

    // matrix for initial conditions
    Vector tCanonicalBasis(m_tFespace->GetVSize());
    tCanonicalBasis = 0.0;
    tCanonicalBasis(0) = 1;
    SparseMatrix *tInit = new SparseMatrix(tCanonicalBasis);

    // compute blocks
    SparseMatrix *tmpBuf;

    // block at row 1, column 1
    tmpBuf = Add(*m_xMass2, *m_xStiff2);
    m_block11 = OuterProduct(*m_tMass, *tmpBuf);
    delete tmpBuf;
    clear(m_xMass2);
    clear(m_xStiff2);

    // block at row 0, column 0
    tmpBuf = Add(*tInit, *m_tStiff);
    m_block00MedIdp = OuterProduct(*tmpBuf, *m_xMass1);
    delete tmpBuf;
    clear(tInit);
    clear(m_tStiff);
    clear(m_xMass1);

    // block at row 0, column 1
    tmpBuf = Transpose(*m_tGrad);
    m_block01MedIdp = OuterProduct(*tmpBuf, *m_xDiv);
    delete tmpBuf;
    clear(m_tGrad);
    clear(m_xDiv);
}

// Assembles medium-dependent system matrices
void heat::LsqXtFemH1Hdiv :: assembleSystemMediumDependent()
{
    // assemble matrices in space

    heat::MediumTensorCoeff medCoeff(m_testCase);

    // xStiff1 form, xV_space
    BilinearForm *xStiff1_form
            = new BilinearForm(m_xFespaces[0]);
    xStiff1_form->AddDomainIntegrator
            (new heat::StiffnessIntegrator(&medCoeff));
    xStiff1_form->Assemble();
    xStiff1_form->Finalize();
    m_xStiff1 = xStiff1_form->LoseMat();
    delete xStiff1_form;

    // xGrad form
    MixedBilinearForm *xGrad_form
            = new MixedBilinearForm(m_xFespaces[0],
                                    m_xFespaces[1]);
    xGrad_form->AddDomainIntegrator
            (new heat::VectorFEGradientIntegrator(&medCoeff));
    xGrad_form->Assemble();
    xGrad_form->Finalize();
    m_xGrad = xGrad_form->LoseMat();
    delete xGrad_form;

    // compute blocks
    SparseMatrix *tmp1, *tmp2, *tmpBuf;

    // block at row 0, column 0
    clear(m_block00);
    tmpBuf = OuterProduct(*m_tMass, *m_xStiff1);
    m_block00 = Add(*m_block00MedIdp, *tmpBuf);
    delete tmpBuf;
    clear(m_xStiff1);

    // block at row 0, column 1
    clear(m_block01);
    tmp2 = Transpose(*m_xGrad);
    tmp1 = OuterProduct(*m_tMass, *tmp2);
    m_block01 = Add(-1, *tmp1, -1, *m_block01MedIdp);
    delete tmp1;
    delete tmp2;
    clear(m_xGrad);

    // block at row 1, column 0
    clear(m_block10);
    m_block10 = Transpose(*m_block01);
}

void heat::LsqXtFemH1Hdiv
:: assembleSourceXDiv(Vector& b, double t) const
{
    LinearForm source_form(m_xFespaces[1]);
    heat::SourceCoeff f(m_testCase);
    f.SetTime(t);
    source_form.AddDomainIntegrator
            (new heat::VectorFEDivergenceLFIntegrator(f,-1));
    source_form.Assemble();
    b = source_form.GetData();
}

// Assembles matrices for L2-projection
void heat::LsqXtFemH1Hdiv :: assembleProjector()
{
    // assemble matrices in time

    // tMass form
    BilinearForm *tMass_form = new BilinearForm(m_tFespace);
    tMass_form->AddDomainIntegrator(new MassIntegrator);
    tMass_form->Assemble();
    tMass_form->Finalize();
    m_tMass = tMass_form->LoseMat();
    delete tMass_form;

    // assemble matrices in space

    // xMass1 form, xV_space
    BilinearForm *xMass1_form
            = new BilinearForm(m_xFespaces[0]);
    xMass1_form->AddDomainIntegrator(new MassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_xMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xMass2 form, xR_space
    BilinearForm *xMass2_form
            = new BilinearForm(m_xFespaces[1]);
    xMass2_form->AddDomainIntegrator
            (new VectorFEMassIntegrator);
    xMass2_form->Assemble();
    xMass2_form->Finalize();
    m_xMass2 = xMass2_form->LoseMat();
    delete xMass2_form;

    //compute blocks
    m_block00 = OuterProduct(*m_tMass, *m_xMass1);
    m_block11 = OuterProduct(*m_tMass, *m_xMass2);
}

void heat::LsqXtFemH1Hdiv
:: assembleXProjectionRhs(Vector& b1, Vector& b2,
                            double t) const
{
    LinearForm uSol_form(m_xFespaces[0]);
    heat::ExactTemperatureCoeff uCoeff(m_testCase);
    uCoeff.SetTime(t);
    uSol_form.AddDomainIntegrator
            (new DomainLFIntegrator(uCoeff));
    uSol_form.Assemble();
    b1 = uSol_form.GetData();

    LinearForm qSol_form(m_xFespaces[1]);
    heat::ExactFluxCoeff qCoeff(m_testCase);
    qCoeff.SetTime(t);
    qSol_form.AddDomainIntegrator
            (new VectorFEDomainLFIntegrator(qCoeff));
    qSol_form.Assemble();
    b2 = qSol_form.GetData();
}

// End of file
