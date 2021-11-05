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
    int order = 2*fe.GetOrder()+1;
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
    elmat = 0.0;

    // set integration rule
    int order = testFeFine.GetOrder() + trialFeCoarse.GetOrder();
    const IntegrationRule *irFine
            = &IntRules.Get(testFeFine.GetGeomType(), order);

    // evaluate integral on the fine element
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
     
// End of file
