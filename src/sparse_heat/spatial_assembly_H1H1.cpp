#include "spatial_assembly.hpp"

using namespace mfem;


// Vector Mass Integrator
void sparseHeat::SpatialVectorMassIntegrator
:: assembleElementMatrix (const FiniteElement & fe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim = fe.GetDim();
    int ndofs = fe.GetDof();

    Vector shape(ndofs);
    DenseMatrix partialElmat(ndofs, ndofs);
    elmat.SetSize(dim*ndofs, dim*ndofs);

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
        MultVVt(shape, partialElmat);
        partialElmat *= weight;

        for (int k = 0; k < dim; k++) {
            elmat.AddMatrix(partialElmat, ndofs*k, ndofs*k);
        }
    }
}

void sparseHeat::SpatialVectorMassIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    Vector trialShapeCoarse(trialNdofsCoarse);
    Vector testShapeFine(testNdofsFine);
    Vector coords(dim);

    DenseMatrix partialElmat(testNdofsFine, trialNdofsCoarse);
    elmat.SetSize(dim*testNdofsFine, dim*trialNdofsCoarse);

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
        testFeFine.CalcShape(ipFine, testShapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialFeCoarse.CalcShape(ipCoarse, trialShapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();

        MultVWt(testShapeFine, trialShapeCoarse, partialElmat);
        partialElmat *= weight;

        for (int k = 0; k < dim; k++) {
            elmat.AddMatrix(partialElmat,
                            testNdofsFine*k, trialNdofsCoarse*k);
        }
    }
}


// Vector Stiffness Integrator
void sparseHeat::SpatialVectorStiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & fe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim  = fe.GetDim();
    int ndofs = fe.GetDof();

    Vector divshape(dim*ndofs);
    DenseMatrix dshape(ndofs, dim);

    elmat.SetSize(dim*ndofs, dim*ndofs);

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

        fe.CalcPhysDShape (elTrans, dshape);
        dshape.GradToDiv (divshape);

        double weight = ip.weight*elTrans.Weight();
        AddMult_a_VVt(weight, divshape, elmat);
    }
}

void sparseHeat::SpatialVectorStiffnessIntegrator
:: assembleElementMatrix (const FiniteElement & testFeFine,
                          ElementTransformation & testElTransFine,
                          const FiniteElement & trialFeCoarse,
                          ElementTransformation & trialElTransCoarse,
                          DenseMatrix & elmat)
{
    int dim = trialFeCoarse.GetDim();
    int testNdofsFine = testFeFine.GetDof();
    int trialNdofsCoarse = trialFeCoarse.GetDof();

    Vector trialDivshapeCoarse(dim*trialNdofsCoarse);
    DenseMatrix trialDshapeCoarse(trialNdofsCoarse, dim);

    Vector testDivshapeFine(dim*testNdofsFine);
    DenseMatrix testDshapeFine(testNdofsFine, dim);

    Vector coords(dim);
    elmat.SetSize(dim*testNdofsFine, dim*trialNdofsCoarse);

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
        testFeFine.CalcPhysDShape (testElTransFine, testDshapeFine);
        testDshapeFine.GradToDiv (testDivshapeFine);

        IntegrationPoint ipCoarse;
        testElTransFine.Transform(ipFine, coords);
        trialElTransCoarse.TransformBack(coords, ipCoarse);
        trialFeCoarse.CalcPhysDShape (trialElTransCoarse, trialDshapeCoarse);
        trialDshapeCoarse.GradToDiv (trialDivshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_VWt(weight, testDivshapeFine, trialDivshapeCoarse, elmat);
    }
}


// Vector Gradient Integrator
void sparseHeat::SpatialVectorGradientIntegrator
:: assembleElementMatrix (const FiniteElement & testFe,
                          const FiniteElement & trialFe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim = testFe.GetDim();
    int testNdofs = testFe.GetDof();
    int trialNdofs = trialFe.GetDof();

    DenseMatrix trialDshape(trialNdofs, dim);
    Vector testShape(testNdofs);
    elmat.SetSize(dim*testNdofs, trialNdofs);

    // set integration rule
    int order = testFe.GetOrder()+trialFe.GetOrder()+2;
    const IntegrationRule *ir
            = &IntRules.Get(testFe.GetGeomType(), order);

    // evaluate integral
    DenseMatrix tmpMat(dim*testNdofs, trialNdofs);
    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        elTrans.SetIntPoint(&ip);

        trialFe.CalcPhysDShape(elTrans, trialDshape);
        testFe.CalcShape(ip, testShape);

        double weight = ip.weight*elTrans.Weight();

        if (m_matrixCoeff) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trialNdofs, dim);
            m_matrixCoeff->Eval(matM, elTrans, ip);
            MultABt(trialDshape, matM, buf);
            assembleBlock(dim, testNdofs, trialNdofs,
                          testShape, buf, tmpMat);
        }
        else {
            assembleBlock(dim, testNdofs, trialNdofs,
                          testShape, trialDshape, tmpMat);
        }

        if (m_scalarCoeff) { // scalar coefficient
            double val = m_scalarCoeff->Eval(elTrans, ip);
            weight *= val;
        }
        elmat.Add(weight, tmpMat);
    }
}

void sparseHeat::SpatialVectorGradientIntegrator
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
    DenseMatrix trialDshapeCoarse(trialNdofsCoarse, dim);
    Vector coords(dim);

    elmat.SetSize(dim*testNdofsFine, trialNdofsCoarse);

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder() + 2;
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    DenseMatrix tmpMat(dim*testNdofsFine, trialNdofsCoarse);
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
        trialFeCoarse.CalcPhysDShape(trialElTransCoarse,
                                     trialDshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();

        if (m_matrixCoeff) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trialNdofsCoarse, dim);
            m_matrixCoeff->Eval(matM, testElTransFine, ipFine);
            MultABt(trialDshapeCoarse, matM, buf);
            assembleBlock(dim, testNdofsFine, trialNdofsCoarse,
                          testShapeFine, buf, tmpMat);
        }
        else {
            assembleBlock(dim, testNdofsFine, trialNdofsCoarse,
                          testShapeFine, trialDshapeCoarse, tmpMat);
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

void sparseHeat::SpatialVectorGradientIntegrator
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
    DenseMatrix trialDshapeFine(trialNdofsFine, dim);
    Vector coords(dim);

    elmat.SetSize(dim*testNdofsCoarse, trialNdofsFine);

    // set integration rule
    int order = testFeCoarse.GetOrder() + trialFeFine.GetOrder() + 2;
    const IntegrationRule *irFine
            = &IntRules.Get(trialFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
    DenseMatrix tmpMat(dim*testNdofsCoarse, trialNdofsFine);
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
        testFeCoarse.CalcShape(ipCoarse, testShapeCoarse);

        double weight = ipFine.weight*trialElTransFine.Weight();

        if (m_matrixCoeff) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trialNdofsFine, dim);
            m_matrixCoeff->Eval(matM, trialElTransFine, ipFine);
            MultABt(trialDshapeFine, matM, buf);
            assembleBlock(dim, testNdofsCoarse, trialNdofsFine,
                          testShapeCoarse, buf, tmpMat);
        }
        else {
            assembleBlock(dim, testNdofsCoarse, trialNdofsFine,
                          testShapeCoarse, trialDshapeFine, tmpMat);
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

void sparseHeat::SpatialVectorGradientIntegrator
:: assembleBlock(const int dim,
                 const int testNdofs,
                 const int trialNdofs,
                 const Vector& shape,
                 const DenseMatrix& dshape,
                 DenseMatrix& outMat) const
{
    for (int k=0; k<dim; k++)
        for (int i=0; i<testNdofs; i++)
            for (int j=0; j<trialNdofs; j++)
            {
                outMat(i + k*testNdofs, j)
                        = shape(i)*dshape(j, k);
            }
}


// Vector Divergence Integrator
void sparseHeat::SpatialVectorDivergenceIntegrator
:: assembleElementMatrix (const FiniteElement & testFe,
                          const FiniteElement & trialFe,
                          ElementTransformation & elTrans,
                          DenseMatrix & elmat)
{
    int dim = trialFe.GetDim();
    int testNdofs = testFe.GetDof();
    int trialNdofs = trialFe.GetDof();

    Vector testShape(testNdofs);
    Vector trialDivshape(dim*trialNdofs);
    DenseMatrix trialDshape(trialNdofs, dim);
    DenseMatrix trialGshape(trialNdofs, dim);
    DenseMatrix trialJadj(dim);

    elmat.SetSize(testNdofs, dim*trialNdofs);

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
        trialFe.CalcPhysDShape (elTrans, trialDshape);
        trialDshape.GradToDiv (trialDivshape);

        double weight = ip.weight * elTrans.Weight();
        AddMult_a_VWt(weight, testShape, trialDivshape, elmat);
    }
}

void sparseHeat::SpatialVectorDivergenceIntegrator
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
    Vector trialDivshapeCoarse(dim*trialNdofsCoarse);
    DenseMatrix trialDshapeCoarse(trialNdofsCoarse, dim);
    Vector coords(dim);

    elmat.SetSize(testNdofsFine, dim*trialNdofsCoarse);

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
        trialFeCoarse.CalcPhysDShape (trialElTransCoarse,
                                      trialDshapeCoarse);
        trialDshapeCoarse.GradToDiv (trialDivshapeCoarse);

        double weight = ipFine.weight*testElTransFine.Weight();
        AddMult_a_VWt(weight, testShapeFine, trialDivshapeCoarse, elmat);
    }
}

void sparseHeat::SpatialVectorDivergenceIntegrator
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
    Vector trialDivshapeFine(dim*trialNdofsFine);
    DenseMatrix trialDshapeFine(trialNdofsFine, dim);
    Vector coords(dim);

    elmat.SetSize(testNdofsCoarse, dim*trialNdofsFine);

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
        trialFeFine.CalcPhysDShape (trialElTransFine,
                                    trialDshapeFine);
        trialDshapeFine.GradToDiv (trialDivshapeFine);

        IntegrationPoint ipCoarse;
        trialElTransFine.Transform(ipFine, coords);
        testElTransCoarse.TransformBack(coords, ipCoarse);
        testFeCoarse.CalcShape(ipCoarse, testShapeCoarse);

        double weight = ipFine.weight * trialElTransFine.Weight();
        AddMult_a_VWt(weight, testShapeCoarse, trialDivshapeFine, elmat);
    }
}
