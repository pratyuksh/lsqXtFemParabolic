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

    Vector testDivshapeFine(testNdofsFine);
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
