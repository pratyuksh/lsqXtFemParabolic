#include "assembly.hpp"
#include <assert.h>


// Mass Integrator
void sparseHeat::MassIntegrator
:: assembleElementMatrix (const FiniteElement & fe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int ndofs = fe.GetDof();

    Vector shape(ndofs);
    elmat.SetSize(ndofs, ndofs);

    // set integration rule
    int order = 2*fe.GetOrder();
    const IntegrationRule *ir
            = &IntRules.Get(fe.GetGeomType(), order);

    // evaluate integral
    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);
        fe.CalcShape(ip, shape);

        double weight = ip.weight*elTrans.Weight();
        AddMult_a_VVt(weight, shape, elmat);
    }
}

void sparseHeat::MassIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    Vector shapeCoarse(trialNdofsCoarse);
    Vector shapeFine(testNdofsFine);
    Vector coords(dim);
    elmat.SetSize(testNdofsFine, trialNdofsCoarse);

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder();
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        testElTransFine.SetIntPoint(&ipFine);
        testFeFine.CalcShape(ipFine, shapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialFeCoarse.CalcShape(ipCoarse, shapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_VWt(weight, shapeFine, shapeCoarse, elmat);
    }
}


// Stiffness Integrator
void sparseHeat::StiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & fe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim = fe.GetDim();
    int ndofs = fe.GetDof();

    DenseMatrix dshape(ndofs, dim);
    elmat.SetSize(ndofs, ndofs);

    // set integration rule
    int order = 2*fe.GetOrder()+2;
    const IntegrationRule *ir
            = &IntRules.Get(fe.GetGeomType(), order);

    // evaluate integral
    elmat = 0.0;
    DenseMatrix tmpMat(ndofs);
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);
        fe.CalcPhysDShape(elTrans, dshape);

        double weight = ip.weight*elTrans.Weight();

        if (m_matrixCoeff) {
            DenseMatrix matM(dim);
            DenseMatrix buf(ndofs, dim);
            m_matrixCoeff->Eval(matM, elTrans, ip);
            MultABt(dshape, matM, buf);
            MultAAt(buf, tmpMat);
        }
        else {
            MultAAt(dshape, tmpMat);
        }

        if (m_scalarCoeff) {
            double val = m_scalarCoeff->Eval(elTrans, ip);
            weight *= val*val;
        }

        elmat.Add(weight, tmpMat);
    }
}

void sparseHeat::StiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    DenseMatrix dshapeCoarse(trialNdofsCoarse, dim);
    DenseMatrix dshapeFine(testNdofsFine, dim);
    Vector coords(dim);

    elmat.SetSize(testNdofsFine, trialNdofsCoarse);
    elmat = 0.0;

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder() + 2;
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    DenseMatrix tmpMat(testNdofsFine, trialNdofsCoarse);
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        testElTransFine.SetIntPoint(&ipFine);
        testFeFine.CalcPhysDShape(testElTransFine, dshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcPhysDShape(trialElTransCoarse, dshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();

        if (m_matrixCoeff) {
            DenseMatrix matM(dim);
            DenseMatrix bufFine(testNdofsFine, dim);
            DenseMatrix bufCoarse(trialNdofsCoarse, dim);
            m_matrixCoeff->Eval(matM, testElTransFine, ipFine);
            MultABt(dshapeFine, matM, bufFine);
            MultABt(dshapeCoarse, matM, bufCoarse);
            MultABt(bufFine, bufCoarse, tmpMat);
        }
        else {
            MultABt(dshapeFine, dshapeCoarse, tmpMat);
        }

        if (m_scalarCoeff) {
            double val = m_scalarCoeff->Eval(testElTransFine, ipFine);
            weight *= val*val;
        }

        elmat.Add(weight, tmpMat);
    }
}


// VectorFE Mass Integrator
void sparseHeat::VectorFEMassIntegrator
:: assembleElementMatrix (const FiniteElement & fe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim = fe.GetDim();
    int ndofs = fe.GetDof();

    DenseMatrix vshape(ndofs, dim);
    elmat.SetSize(ndofs, ndofs);

    // set integration rule
    int order = 2*fe.GetOrder();
    const IntegrationRule *ir
            = &IntRules.Get(fe.GetGeomType(), order);

    // evaluate integral
    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);
        fe.CalcVShape(elTrans, vshape);

        double weight = ip.weight*elTrans.Weight();
        AddMult_a_AAt(weight, vshape, elmat);
    }
}

void sparseHeat::VectorFEMassIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    DenseMatrix vshapeCoarse(trialNdofsCoarse, dim);
    DenseMatrix vshapeFine(testNdofsFine, dim);
    Vector coords(dim);
    elmat.SetSize(testNdofsFine, trialNdofsCoarse);

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder();
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        testElTransFine.SetIntPoint(&ipFine);
        testFeFine.CalcVShape(testElTransFine, vshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcVShape(trialElTransCoarse, vshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_ABt(weight, vshapeFine, vshapeCoarse, elmat);
    }
}


// VectorFE Stiffness Integrator
void sparseHeat::VectorFEStiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & fe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int ndofs = fe.GetDof();

    Vector divshape(ndofs);
    elmat.SetSize(ndofs, ndofs);

    // set integration rule
    int order = 2*fe.GetOrder()-2;
    const IntegrationRule *ir
            = &IntRules.Get(fe.GetGeomType(), order);

    // evaluate integral
    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);
        fe.CalcPhysDivShape(elTrans, divshape);

        double weight = ip.weight*elTrans.Weight();
        AddMult_a_VVt(weight, divshape, elmat);
    }
}

void sparseHeat::VectorFEStiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    Vector divshapeCoarse(trialNdofsCoarse);
    Vector divshapeFine(testNdofsFine);
    Vector coords(dim);
    elmat.SetSize(testNdofsFine, trialNdofsCoarse);

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder() - 2;
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        testElTransFine.SetIntPoint(&ipFine);
        testFeFine.CalcPhysDivShape(testElTransFine, divshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcPhysDivShape(trialElTransCoarse,
                                       divshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_VWt(weight, divshapeFine, divshapeCoarse, elmat);
    }
}

// End of file
