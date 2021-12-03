#include "spatial_assembly.hpp"

using namespace mfem;


// VectorFE Mass Integrator
void sparseHeat::SpatialVectorFEMassIntegrator
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

void sparseHeat::SpatialVectorFEMassIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    DenseMatrix trialVshapeCoarse(trialNdofsCoarse, dim);
    DenseMatrix testVshapeFine(testNdofsFine, dim);
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
        testFeFine.CalcVShape(testElTransFine, testVshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcVShape(trialElTransCoarse, trialVshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_ABt(weight, testVshapeFine, trialVshapeCoarse, elmat);
    }
}


// VectorFE Stiffness Integrator
void sparseHeat::SpatialVectorFEStiffnessIntegrator
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

void sparseHeat::SpatialVectorFEStiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    Vector trialDivshapeCoarse(trialNdofsCoarse);
    Vector testDivshapeFine(testNdofsFine);
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
        testFeFine.CalcPhysDivShape(testElTransFine, testDivshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcPhysDivShape(trialElTransCoarse,
                                       trialDivshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_VWt(weight, testDivshapeFine, trialDivshapeCoarse, elmat);
    }
}


// VectorFE Gradient Integrator
void sparseHeat::SpatialVectorFEGradientIntegrator
:: assembleElementMatrix (const FiniteElement & testFe,
                          const FiniteElement & trialFe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim = testFe.GetDim();
    int testNdofs = testFe.GetDof();
    int trialNdofs = trialFe.GetDof();

    DenseMatrix trialDshape(trialNdofs, dim);
    DenseMatrix testVshape(testNdofs, dim);
    elmat.SetSize(testNdofs, trialNdofs);

    // set integration rule
    int order = testFe.GetOrder()+trialFe.GetOrder()+2;
    const IntegrationRule *ir
            = &IntRules.Get(testFe.GetGeomType(), order);

    // evaluate integral
    DenseMatrix tmpMat(testNdofs, trialNdofs);
    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);

        trialFe.CalcPhysDShape(elTrans, trialDshape);
        testFe.CalcVShape(elTrans, testVshape);

        double weight = ip.weight*elTrans.Weight();

        if (m_matrixCoeff) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trialNdofs, dim);
            m_matrixCoeff->Eval(matM, elTrans, ip);
            MultABt(trialDshape, matM, buf);
            MultABt(testVshape, buf, tmpMat);
        }
        else {
            MultABt(testVshape, trialDshape, tmpMat);
        }

        if (m_scalarCoeff) { // scalar coefficient
            double val = m_scalarCoeff->Eval(elTrans, ip);
            weight *= val;
        }
        elmat.Add(weight, tmpMat);
    }
}

void sparseHeat::SpatialVectorFEGradientIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    DenseMatrix testVshapeFine(testNdofsFine, dim);
    DenseMatrix trialDshapeCoarse(trialNdofsCoarse, dim);
    Vector coords(dim);
    elmat.SetSize(testNdofsFine, trialNdofsCoarse);

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
        testFeFine.CalcVShape(testElTransFine, testVshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcPhysDShape(trialElTransCoarse,
                                     trialDshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();

        if (m_matrixCoeff) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trialNdofsCoarse, dim);
            m_matrixCoeff->Eval(matM, testElTransFine, ipFine);
            MultABt(trialDshapeCoarse, matM, buf);
            MultABt(testVshapeFine, buf, tmpMat);
        }
        else {
            MultABt(testVshapeFine, trialDshapeCoarse, tmpMat);
        }

        if (m_scalarCoeff) { // scalar coefficient
            double val
                    = m_scalarCoeff->Eval(testElTransFine,
                                          ipFine);
            weight *= val;
        }
        elmat.Add(weight, tmpMat);
    }
}

void sparseHeat::SpatialVectorFEGradientIntegrator
:: assembleElementMatrix2 (const FiniteElement & testFeCoarse,
                           ElementTransformation & testElTransCoarse,
                           const FiniteElement & trialFeFine,
                           ElementTransformation & trialElTransFine,
                           DenseMatrix & elmat)
{
    int dim = trialFeFine.GetDim();
    int testNdofsCoarse = testFeCoarse.GetDof();
    int trialNdofsFine = trialFeFine.GetDof();

    DenseMatrix testVshapeCoarse(testNdofsCoarse, dim);
    DenseMatrix trialDshapeFine(trialNdofsFine, dim);
    Vector coords(dim);
    elmat.SetSize(testNdofsCoarse, trialNdofsFine);

    // set integration rule
    int order = testFeCoarse.GetOrder() + trialFeFine.GetOrder() + 2;
    const IntegrationRule *irFine
            = &IntRules.Get(trialFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    DenseMatrix tmpMat(testNdofsCoarse, trialNdofsFine);
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        trialElTransFine.SetIntPoint(&ipFine);
        trialFeFine.CalcPhysDShape(trialElTransFine,
                                   trialDshapeFine);

        IntegrationPoint ipCoarse;
        trialElTransFine.Transform(ipFine, coords);
        testElTransCoarse.TransformBack(coords, ipCoarse);
        testElTransCoarse.SetIntPoint(&ipCoarse);
        testFeCoarse.CalcVShape(testElTransCoarse,
                                testVshapeCoarse);

        double weight = ipFine.weight*trialElTransFine.Weight();

        if (m_matrixCoeff) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trialNdofsFine, dim);
            m_matrixCoeff->Eval(matM, trialElTransFine, ipFine);
            MultABt(trialDshapeFine, matM, buf);
            MultABt(testVshapeCoarse, buf, tmpMat);
        }
        else {
            MultABt(testVshapeCoarse, trialDshapeFine, tmpMat);
        }

        if (m_scalarCoeff) { // scalar coefficient
            double val
                    = m_scalarCoeff->Eval(trialElTransFine,
                                          ipFine);
            weight *= val;
        }
        elmat.Add(weight, tmpMat);
    }
}


// VectorFE Divergence Integrator
void sparseHeat::SpatialVectorFEDivergenceIntegrator
:: assembleElementMatrix (const FiniteElement & testFe,
                          const FiniteElement & trialFe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int testNdofs = testFe.GetDof();
    int trialNdofs = trialFe.GetDof();

    Vector testShape(testNdofs);
    Vector trialDivshape(trialNdofs);
    elmat.SetSize(testNdofs, trialNdofs);

    // set integration rule
    int order = testFe.GetOrder()+trialFe.GetOrder()-1;
    const IntegrationRule *ir
            = &IntRules.Get(testFe.GetGeomType(), order);

    // evaluate integral
    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);

        testFe.CalcShape(ip, testShape);
        trialFe.CalcDivShape(ip, trialDivshape);

        double weight = ip.weight;
        AddMult_a_VWt(weight, testShape, trialDivshape, elmat);
    }
}

void sparseHeat::SpatialVectorFEDivergenceIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    Vector testShapeFine(testNdofsFine);
    Vector trialDivshapeCoarse(trialNdofsCoarse);
    Vector coords(dim);
    elmat.SetSize(testNdofsFine, trialNdofsCoarse);

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder() - 1;
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        testElTransFine.SetIntPoint(&ipFine);
        testFeFine.CalcShape(ipFine, testShapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialElTransCoarse.SetIntPoint(&ipCoarse);
        trialFeCoarse.CalcPhysDivShape(trialElTransCoarse,
                                       trialDivshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_VWt(weight, testShapeFine, trialDivshapeCoarse, elmat);
    }
}

void sparseHeat::SpatialVectorFEDivergenceIntegrator
:: assembleElementMatrix2 (const FiniteElement & testFeCoarse,
                           ElementTransformation & testElTransCoarse,
                           const FiniteElement & trialFeFine,
                           ElementTransformation & trialElTransFine,
                           DenseMatrix & elmat)
{
    int dim = trialFeFine.GetDim();
    int testNdofsCoarse = testFeCoarse.GetDof();
    int trialNdofsFine = trialFeFine.GetDof();

    Vector testShapeCoarse(testNdofsCoarse);
    Vector testDivshapeFine(trialNdofsFine);
    Vector coords(dim);
    elmat.SetSize(testNdofsCoarse, trialNdofsFine);

    // set integration rule
    int order = testFeCoarse.GetOrder() + trialFeFine.GetOrder() - 1;
    const IntegrationRule *irFine
            = &IntRules.Get(trialFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    elmat = 0.0;
    for (int i = 0; i < irFine->GetNPoints(); i++)
    {
        const IntegrationPoint &ipFine = irFine->IntPoint(i);
        trialElTransFine.SetIntPoint(&ipFine);
        trialFeFine.CalcDivShape(ipFine, testDivshapeFine);

        IntegrationPoint ipCoarse;
        trialElTransFine.Transform(ipFine, coords);
        testElTransCoarse.TransformBack(coords, ipCoarse);
        testFeCoarse.CalcShape(ipCoarse, testShapeCoarse);

        double weight = ipFine.weight;
        AddMult_a_VWt(weight, testShapeCoarse, testDivshapeFine, elmat);
    }
}

// End of file
