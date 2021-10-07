#include "functionals.hpp"
#include <assert.h>


// Quantity of interest LF Integrator
void heat::QoILFIntegrator
:: AssembleRHSElementVect (const FiniteElement &fe,
                           ElementTransformation &Trans,
                           Vector &elvect)
{
    int dim = fe.GetDim();
    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(ndofs, dim);
#else
    dshape.SetSize(ndofs, dim);
#endif
    elvect.SetSize(ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = fe.GetOrder()+4;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    // compute x^(m_pow) element-wise
    DenseMatrix xPow(dim, ir->GetNPoints());
    Trans.Transform(*ir, xPow);
    for (int j=0; j < ir->GetNPoints(); j++) {
        for (int i=0; i<dim; i++) {
            xPow(i,j) = std::pow(xPow(i,j), m_pow);
        }
    }

    elvect = 0.0;
    Vector v(ndofs);
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        Vector xPow_;
        xPow.GetColumnReference(i, xPow_);

        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);
        fe.CalcPhysDShape(Trans, dshape);

        double w = ip.weight*Trans.Weight();
        dshape.Mult(xPow_, v);
        elvect.Add(w, v);
    }
}

// Quantity of interest functional
void heat::QoIFunctional
:: init () const
{
    m_bObs.SetSize(m_fes->GetTrueVSize());

    LinearForm lf(m_fes);
    lf.AddDomainIntegrator(new heat::QoILFIntegrator(m_pow));
    lf.Assemble();
    m_bObs = lf.GetData();
    //std::cout << "\n\n";
    //m_bObs.Print();
}

// Evalautes the QoI for a function passed as an MFEM GridFunction
double heat::QoIFunctional
:: operator()(const GridFunction &u) const
{
    //std::cout << "\n\n";
    //u.Print();
    return std::move(u*m_bObs);
}

// Evalautes the QoI for a function passed as a vector of coefficients
double heat::QoIFunctional
:: operator()(const Vector &U) const
{
    //std::cout << "\n\n";
    //u.Print();
    return std::move(U*m_bObs);
}

// Quantity of interest xt functional
void heat::QoIXtFunctional
:: init () const
{
    int tDim = m_tFes->GetTrueVSize();
    int xDim = m_xFes->GetTrueVSize();

    // x functional
    Vector bx(xDim);
    LinearForm lf(m_xFes);
    lf.AddDomainIntegrator(new heat::QoILFIntegrator(m_pow));
    lf.Assemble();
    bx = lf.GetData();

    int Nt = m_tFes->GetNE();
    m_bObs.SetSize(tDim*xDim);
    m_bObs = 0.;

    Array<int> vdofs;
    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector elvec, shape;
    for (int n=0; n<Nt; n++)
    {
        m_tFes->GetElementVDofs (n, vdofs);
        tTrans = m_tFes->GetElementTransformation(n);

        // compute elvec
        tFe = m_tFes->GetFE(n);
        int tNdofs = tFe->GetDof();
        shape.SetSize(tNdofs);
        elvec.SetSize(tNdofs*xDim);

        int order = 2*tFe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(tFe->GetGeomType(), order);

        elvec = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            tTrans->SetIntPoint(&ip);
            tFe->CalcShape(ip, shape);

            // elvec += w*shape*bx
            double w = ip.weight*tTrans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*shape(j);
                for (int k=0; k<xDim; k++) {
                    elvec(k + j*xDim) += coeff*bx(k);
                }
            }
        }
        add_vector(vdofs, elvec, m_bObs);
    }
}

void heat::QoIXtFunctional
:: add_vector(const Array<int> &vdofs,
              const Vector& elvec,
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

// Evaluate the QoI
double heat::QoIXtFunctional
:: operator()(const Vector &U) const
{
    //std::cout << "\n\n";
    //u.Print();
    return std::move(U*m_bObs);
}

// End of file
