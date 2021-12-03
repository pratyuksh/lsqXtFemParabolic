#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"

#include <fstream>

using namespace mfem;


// Constructor
heat::LsqXtFEM :: LsqXtFEM (const nlohmann::json& config,
                            std::shared_ptr<heat::TestCases>& testCase)
    : m_config (config),
      m_testCase(testCase)
{
    m_deg = config["deg"];
}

// Destructor
heat::LsqXtFEM :: ~LsqXtFEM()
{
    reset();
}

// Releases all the memory allocated by this class
void heat::LsqXtFEM :: reset()
{
    m_firstPass = true;

    if (m_tFec) { delete m_tFec; }
    if (m_tFespace) { delete m_tFespace; }

    for (int i=0; i<m_xFecs.Size(); i++)
    {
        if (m_xFecs[i]) { delete m_xFecs[i]; }
        if (m_xFespaces[i]) { delete m_xFespaces[i]; }
    }

    if (m_heatOp) { delete m_heatOp; }
    if (m_heatMat) { delete m_heatMat; }
    if (m_heatProjMat) { delete m_heatProjMat; }

    if (m_tMass) { clear(m_tMass); }
    if (m_tStiff) { clear(m_tStiff); }
    if (m_tGrad) { clear(m_tGrad); }

    if (m_xMass1) { clear(m_xMass1); }
    if (m_xMass2) { clear(m_xMass2); }
    if (m_xStiff1) { clear(m_xStiff1); }
    if (m_xStiff2) { clear(m_xStiff2); }
    if (m_xGrad) { clear(m_xGrad); }
    if (m_xDiv) { clear(m_xDiv); }

    if (m_block00) { clear(m_block00); }
    if (m_block01) { clear(m_block01); }
    if (m_block10) { clear(m_block10); }
    if (m_block11) { clear(m_block11); }

    if (m_block00MedIdp) { clear(m_block00MedIdp); }
    if (m_block01MedIdp) { clear(m_block01MedIdp); }
}

// Assembles the blocks of the linear system matrix
void heat::LsqXtFEM :: assembleSystem()
{   
    if (m_firstPass) {
        assembleSystemMediumIndependent();
        m_firstPass = false;
    }
    assembleSystemMediumDependent();
}

// Builds the linear system matrix from the blocks as an MFEM operator
void heat::LsqXtFEM :: buildSystemOp()
{
    if (m_heatOp) {
        delete m_heatOp;
        m_heatOp = nullptr;
    }
    m_heatOp = new BlockOperator(m_block_offsets);
    m_heatOp->SetBlock(0,0, m_block00);
    m_heatOp->SetBlock(0,1, m_block01);
    m_heatOp->SetBlock(1,0, m_block10);
    m_heatOp->SetBlock(1,1, m_block11);
}

// Builds the linear system matrix from the blocks as a monolithic matrix
void heat::LsqXtFEM :: buildSystemMatrix()
{
    BlockMatrix *heatBlockMat
            = new BlockMatrix(m_block_offsets);
    heatBlockMat->SetBlock(0,0, m_block00);
    heatBlockMat->SetBlock(0,1, m_block01);
    heatBlockMat->SetBlock(1,0, m_block10);
    heatBlockMat->SetBlock(1,1, m_block11);

    if (m_heatMat) { clear(m_heatMat); }
    m_heatMat = heatBlockMat->CreateMonolithic();
    delete heatBlockMat;

    applyBCs(*m_heatMat);

    if (!m_heatMat->ColumnsAreSorted()) {
        m_heatMat->SortColumnIndices();
    }
}

// Builds system matrix from heatOp
// stores only the upper triangle
void heat::LsqXtFEM :: buildSystemMatrixUpperTriangle()
{
    BlockMatrix *heatBlockMat
            = new BlockMatrix(m_block_offsets);
    heatBlockMat->SetBlock(0,0, m_block00);
    heatBlockMat->SetBlock(0,1, m_block01);
    heatBlockMat->SetBlock(1,0, m_block10);
    heatBlockMat->SetBlock(1,1, m_block11);

    SparseMatrix *heatMatBuf = heatBlockMat->CreateMonolithic();
    delete heatBlockMat;

    applyBCs(*heatMatBuf);
    if (m_heatMat) { clear(m_heatMat); }
    m_heatMat = &getUpperTriangle(*heatMatBuf);
    delete heatMatBuf;

    if (!m_heatMat->ColumnsAreSorted()) {
        m_heatMat->SortColumnIndices();
    }
}

// Assembles rhs
void heat::LsqXtFEM :: assembleRhs(BlockVector* B) const
{    
    Vector &B0 = B->GetBlock(0);
    assembleICs(B0);
    assembleSourceTGradXProjection(B0);

    Vector &B1 = B->GetBlock(1);
    assembleSourceTProjectionXDiv(B1);

    applyBCs(*B);
}

// Assembles initial conditions
void heat::LsqXtFEM :: assembleICs(Vector& b) const
{
    int xdimV = m_xFespaces[0]->GetTrueVSize();
    Vector bx(xdimV);
    assembleICsXProjection(bx);

    Array<int> vdofs;
    const FiniteElement *fe = nullptr;
    Vector elvec, shape;

    int n=0; // only for element (0, dt)
    {
        m_tFespace->GetElementVDofs (n, vdofs);
        assert(vdofs[0] == 0);

        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec.SetSize(tNdofs*xdimV);

        IntegrationPoint ip;
        ip.Set1w(0.0, 1.0);
        fe->CalcShape(ip, shape);

        for (int j=0; j<tNdofs; j++) {
            for (int k=0; k<xdimV; k++) {
                elvec(k + j*xdimV) = shape(j)*bx(k);
            }
        }
        addVector(vdofs, elvec, b);
    }
}

void heat::LsqXtFEM
:: assembleICsXProjection(Vector& b) const
{
    LinearForm ics_form(m_xFespaces[0]);
    heat::InitialTemperatureCoeff f(m_testCase);
    ics_form.AddDomainIntegrator
            (new DomainLFIntegrator(f));
    ics_form.Assemble();
    b = ics_form.GetData();
}

// Assembles source
void heat::LsqXtFEM
:: assembleSourceTGradXProjection(Vector &b) const
{
    int Nt = m_tFespace->GetNE();
    int xdimV = m_xFespaces[0]->GetTrueVSize();
    Vector bx(xdimV);

    Array<int> vdofs;
    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec;
    DenseMatrix dshape;

    for (int n=0; n<Nt; n++)
    {
        m_tFespace->GetElementVDofs (n, vdofs);
        trans = m_tFespace->GetElementTransformation(n);

        // compute elvec
        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        dshape.SetSize(tNdofs, 1); // time dimension is 1d
        elvec.SetSize(tNdofs*xdimV);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcPhysDShape(*trans, dshape);

            Vector t;
            trans->Transform(ip, t);
            assembleSourceXProjection(bx, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*dshape(j,0);
                for (int k=0; k<xdimV; k++) {
                    elvec(k + j*xdimV) += coeff*bx(k);
                }
            }
        }
        addVector(vdofs, elvec, b);
    }
}

void heat::LsqXtFEM
:: assembleSourceXProjection(Vector& b,
                               double t) const
{
    LinearForm source_form(m_xFespaces[0]);
    heat::SourceCoeff f(m_testCase);
    f.SetTime(t);
    source_form.AddDomainIntegrator
            (new DomainLFIntegrator(f, 2, 4));
    source_form.Assemble();
    b = source_form.GetData();
}

void heat::LsqXtFEM
:: assembleSourceTProjectionXDiv(Vector& b) const
{
    int Nt = m_tFespace->GetNE();
    int xdimR = m_xFespaces[1]->GetTrueVSize();
    Vector bx(xdimR);

    Array<int> vdofs;
    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec, shape;

    for (int n=0; n<Nt; n++)
    {
        m_tFespace->GetElementVDofs (n, vdofs);
        trans = m_tFespace->GetElementTransformation(n);

        // compute elvec
        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec.SetSize(tNdofs*xdimR);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            Vector t;
            trans->Transform(ip, t);
            assembleSourceXDiv(bx, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*shape(j);
                for (int k=0; k<xdimR; k++) {
                    elvec(k + j*xdimR) += coeff*bx(k);
                }
            }
        }
        addVector(vdofs, elvec, b);
    }
}

// Adds elvec to b
void heat::LsqXtFEM
:: addVector(const Array<int> &vdofs, const Vector& elvec,
             Vector& b) const
{
    int ndofs = elvec.Size();
    int tNdofs = vdofs.Size();
    int xNdofs = ndofs/tNdofs;

    const bool use_dev
            = vdofs.UseDevice() || elvec.UseDevice();
    auto d_elvec = elvec.Read(use_dev);
    auto d_b = b.ReadWrite(use_dev);

    for(int i=0; i<tNdofs; i++)
    {
        const int j = vdofs[i];
        const int ii = i*xNdofs;
        if (j >= 0) {
            const int jj = j*xNdofs;
            for (int k=0; k<xNdofs; k++)
                d_b[k + jj] += d_elvec[k + ii];
        }
        else {
            const int jj = (-1-j)*xNdofs;
            for (int k=0; k<xNdofs; k++)
                d_b[k + jj] -= d_elvec[k + ii];
        }
    }
}

// Applies BCs to SparseMatrix
void heat::LsqXtFEM :: applyBCs(SparseMatrix& A) const
{
    for (int k=0; k<m_essTdofList.Size(); k++) {
        A.EliminateRowCol(m_essTdofList[k]);
    }
}

// Applies BCs to BlockVector
void heat::LsqXtFEM :: applyBCs(BlockVector& B) const
{
    Vector& B0 = B.GetBlock(0);
    B0.SetSubVector(m_essTdofList, 0.0);
}


// Builds projection matrix
// stores only the upper triangle
void heat::LsqXtFEM :: buildProjectorMatrix()
{
    BlockMatrix *heatProjBlockMat
            = new BlockMatrix(m_block_offsets);
    heatProjBlockMat->SetBlock(0,0, m_block00);
    heatProjBlockMat->SetBlock(1,1, m_block11);

    SparseMatrix *heatProjMatBuf = heatProjBlockMat->CreateMonolithic();
    applyBCs(*heatProjMatBuf);
    m_heatProjMat = &getUpperTriangle(*heatProjMatBuf);
    delete heatProjBlockMat;
    delete heatProjMatBuf;

    if (!m_heatProjMat->ColumnsAreSorted()) {
        m_heatProjMat->SortColumnIndices();
    }
}

// Assembles projection rhs
void heat::LsqXtFEM
:: assembleProjectionRhs(BlockVector* B) const
{
    int Nt = m_tFespace->GetNE();
    int xdimV = m_xFespaces[0]->GetTrueVSize();
    int xdimR = m_xFespaces[1]->GetTrueVSize();
    Vector bx1(xdimV), bx2(xdimR);

    Array<int> vdofs;
    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec1, elvec2, shape;

    for (int n=0; n<Nt; n++)
    {
        m_tFespace->GetElementVDofs (n, vdofs);
        trans = m_tFespace->GetElementTransformation(n);

        // compute elvec
        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec1.SetSize(tNdofs*xdimV);
        elvec2.SetSize(tNdofs*xdimR);

        int order = 2*fe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec1 = 0.0;
        elvec2 = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            Vector t;
            trans->Transform(ip, t);
            assembleXProjectionRhs(bx1, bx2, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*shape(j);
                for (int k=0; k<xdimV; k++) {
                    elvec1(k + j*xdimV) += coeff*bx1(k);
                }
                for (int k=0; k<xdimR; k++) {
                    elvec2(k + j*xdimR) += coeff*bx2(k);
                }
            }
        }
        addVector(vdofs, elvec1, B->GetBlock(0));
        addVector(vdofs, elvec2, B->GetBlock(1));
    }
    applyBCs(*B);
}

// End of file
