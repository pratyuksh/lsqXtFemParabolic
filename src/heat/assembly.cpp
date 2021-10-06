#include "assembly.hpp"
#include <assert.h>


// Gradient Integrator in time, dim = 1
void heat::GradientIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &elmat)
{
    int dim = fe.GetDim();
    assert(dim == 1);

    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector shape(ndofs);
    DenseMatrix dshape(ndofs, dim);
#else
    shape.SetSize(ndofs);
    dshape.SetSize(ndofs, dim);
#endif
    elmat.SetSize(ndofs, ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe.GetOrder()-1;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcShape(ip, shape);
        fe.CalcPhysDShape(Trans, dshape);
        double w = ip.weight*Trans.Weight();

        Vector t;
        Trans.Transform(ip,t);

        for (int l=0; l<ndofs; l++) {
            for (int k=0; k<ndofs; k++) {
                elmat(k,l) += w*shape(k)*dshape(l,0);
            }
        }
    }
}


// Stiffness Integrator
void heat::StiffnessIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &elmat)
{
    int dim = fe.GetDim();
    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(ndofs, dim);
#else
    dshape.SetSize(ndofs, dim);
#endif
    elmat.SetSize(ndofs, ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        //int order = 2*fe.GetOrder()-2;
        int order = 2*fe.GetOrder()+4;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcPhysDShape(Trans, dshape);
        double w = ip.weight*Trans.Weight();

        DenseMatrix tmpMat(ndofs);
        if (M) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(ndofs, dim);
            M->Eval(matM, Trans, ip);
            MultABt(dshape, matM, buf);
            MultAAt(buf, tmpMat);
        }
        else {
            MultAAt(dshape, tmpMat);
        }

        if (q) { // scalar coefficient
            double val = q->Eval(Trans, ip);
            w *= val*val;
        }

        elmat.Add(w,tmpMat);
    }
}


// VectorFE Stiffness Integrator
// Uses, for example, Raviart-Thomas space
void heat::VectorFEStiffnessIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &elmat)
{
    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector divshape(ndofs);
#else
    divshape.SetSize(ndofs);
#endif
    elmat.SetSize(ndofs, ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe.GetOrder()-2;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcPhysDivShape(Trans, divshape);
        double w = ip.weight*Trans.Weight();

        AddMult_a_VVt(w, divshape, elmat);

        /*for (int l=0; l<ndofs; l++) {
            for (int k=0; k<ndofs; k++) {
                elmat(k,l) += w*divshape(k)*divshape(l);
            }
        }*/
    }
}


// VectorFE Gradient Integrator
// Uses, for example, Raviart-Thomas space
void heat::VectorFEGradientIntegrator
:: AssembleElementMatrix2 (const FiniteElement &trial_fe,
                           const FiniteElement &test_fe,
                           ElementTransformation &Trans,
                           DenseMatrix &elmat)
{
    int dim = trial_fe.GetDim();
    int trial_ndofs = trial_fe.GetDof();
    int test_ndofs = test_fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    DenseMatrix gradshape(trial_ndofs, dim);
    DenseMatrix vshape(test_ndofs, dim);
#else
    gradshape.SetSize(trial_ndofs, dim);
    vshape.SetSize(test_ndofs, dim);
#endif
    elmat.SetSize(test_ndofs, trial_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        //int order = trial_fe.GetOrder()+test_fe.GetOrder()-1;
        int order = trial_fe.GetOrder()+test_fe.GetOrder()+4;
        ir = &IntRules.Get(trial_fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        trial_fe.CalcPhysDShape(Trans, gradshape);
        test_fe.CalcVShape(Trans, vshape);
        double w = ip.weight*Trans.Weight();

        DenseMatrix tmpMat(test_ndofs, trial_ndofs);
        if (M) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trial_ndofs, dim);
            M->Eval(matM, Trans, ip);
            MultABt(gradshape, matM, buf);
            MultABt(vshape, buf, tmpMat);
        }
        else {
            MultABt(vshape, gradshape, tmpMat);
        }

        if (q) { // scalar coefficient
            double val = q->Eval(Trans, ip);
            w *= val;
        }
        elmat.Add(w,tmpMat);
    }
}


// VectorFE Divergence LF Integrator
// Uses, for example, Raviart-Thomas space
void heat::VectorFEDivergenceLFIntegrator
:: AssembleRHSElementVect (const FiniteElement &fe,
                           ElementTransformation &Trans,
                           Vector &elvect)
{
    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector divshape(ndofs);
#else
    divshape.SetSize(ndofs);
#endif
    elvect.SetSize(ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe.GetOrder()+2;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcPhysDivShape(Trans, divshape);
        double w = ip.weight*Trans.Weight();
        w *= F.Eval(Trans, ip);

        elvect.Add(w, divshape);

        /*for (int k=0; k<ndofs; k++) {
            elvect(k) += w*divshape(k);
        }*/
    }
    elvect *= m_coeff;
}


// Vector Stiffness Integrator
void heat::VectorStiffnessIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &elmat)
{
    int dim  = fe.GetDim();
    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector divshape(dim*ndofs);
    DenseMatrix dshape(ndofs, dim);
#else
    divshape.SetSize(dim*ndofs);
    dshape.SetSize (ndofs, dim);
#endif
    elmat.SetSize(dim*ndofs, dim*ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe.GetOrder()-2;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcPhysDShape (Trans, dshape);
        dshape.GradToDiv (divshape);

        double w = ip.weight*Trans.Weight();
        AddMult_a_VVt(w, divshape, elmat);
    }
}


// Vector Gradient Integrator
void heat::VectorGradientIntegrator
:: AssembleElementMatrix2 (const FiniteElement &trial_fe,
                           const FiniteElement &test_fe,
                           ElementTransformation &Trans,
                           DenseMatrix &elmat)
{
    int dim = trial_fe.GetDim();
    int trial_ndofs = trial_fe.GetDof();
    int test_ndofs = test_fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(trial_ndofs, dim);
    Vector shape(test_ndofs);
#else
    dshape.SetSize(trial_ndofs, dim);
    shape.SetSize(test_ndofs);
#endif
    elmat.SetSize(test_ndofs*dim, trial_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order
                = trial_fe.GetOrder()+test_fe.GetOrder()-1;
        ir = &IntRules.Get(trial_fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        trial_fe.CalcPhysDShape(Trans, dshape);
        test_fe.CalcShape(ip, shape);
        double w = ip.weight*Trans.Weight();

        DenseMatrix tmpMat(dim*test_ndofs, trial_ndofs);
        if (M) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trial_ndofs, dim);
            M->Eval(matM, Trans, ip);
            MultABt(dshape, matM, buf);
            AssembleBlock(dim, test_ndofs, trial_ndofs,
                          shape, buf, tmpMat);
        }
        else {
            AssembleBlock(dim, test_ndofs, trial_ndofs,
                          shape, dshape, tmpMat);
        }

        if (q) { // scalar coefficient
            double val = q->Eval(Trans, ip);
            w *= val;
        }
        elmat.Add(w,tmpMat);
    }
}

// Vector Gradient Integrator
void heat::VectorGradientIntegrator
:: AssembleBlock(const int dim,
                 const int test_ndofs,
                 const int trial_ndofs,
                 const Vector& shape,
                 const DenseMatrix& dshape,
                 DenseMatrix& outMat) const
{
    for (int k=0; k<dim; k++)
        for (int i=0; i<test_ndofs; i++)
            for (int j=0; j<trial_ndofs; j++)
            {
                outMat(i + k*test_ndofs, j)
                        = shape(i)*dshape(j,k);
            }
}

// Vector Divergence LF Integrator
void heat::VectorDivergenceLFIntegrator
:: AssembleRHSElementVect (const FiniteElement &fe,
                           ElementTransformation &Trans,
                           Vector &elvect)
{
    int dim  = fe.GetDim();
    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector divshape(dim*ndofs);
    DenseMatrix dshape(ndofs, dim);
#else
    divshape.SetSize(dim*ndofs);
    dshape.SetSize (ndofs, dim);
#endif
    elvect.SetSize(dim*ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe.GetOrder()+2;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcPhysDShape (Trans, dshape);
        dshape.GradToDiv (divshape);

        double w = ip.weight*Trans.Weight();
        w *= F.Eval(Trans, ip);
        elvect.Add(w, divshape);
    }
    elvect *= m_coeff;
}

// End of file
