#include "spatial_error_evaluator.hpp"

using namespace mfem;


void sparseHeat::SpatialErrorOfSolution
:: setCoefficients(std::shared_ptr<Coefficient> & uE,
                   std::shared_ptr<VectorCoefficient> & qE,
                   std::shared_ptr<VectorCoefficient> & dudxE,
                   std::shared_ptr<Coefficient> & dudtE,
                   std::shared_ptr<Coefficient> & source)
{
        m_exactTemperatureCoeff = uE;
        m_exactHeatFluxCoeff = qE;

        m_exactSpatialGradientOfTemperatureCoeff = dudxE;
        m_exactTemporalGradientOfTemperatureCoeff = dudtE;

        m_sourceCoeff = source;
}

void sparseHeat::SpatialErrorOfSolution
:: setCurrentTimeForCoefficients(double t) const
{
        m_exactTemperatureCoeff->SetTime(t);
        m_exactHeatFluxCoeff->SetTime(t);

        m_exactSpatialGradientOfTemperatureCoeff->SetTime(t);
        m_exactTemporalGradientOfTemperatureCoeff->SetTime(t);

        m_sourceCoeff->SetTime(t);
}

void sparseHeat::SpatialErrorOfSolution
:: initializeSpatialParentElementIds(int k, Array<int> &elIds) const
{
    assert(elIds.Size() == m_numLevels+1);

    auto hierMeshTrans = m_spatialMeshHierarchy->getTransformations();

    elIds[0] = k;
    elIds[1] = k;

    for (int i=1; i<m_numLevels; i++) {
        elIds[i+1]
                = hierMeshTrans[m_numLevels-1-i]->getParentId(elIds[i]);
    }
}

void sparseHeat::SpatialErrorOfSolution
:: getSpatialParentElementIds(int k, Array<int> &elIds) const
{
    assert(elIds.Size() == m_numLevels+1);

    auto hierMeshTrans = m_spatialMeshHierarchy->getTransformations();
    Array<int> oldElIds(elIds);
//    oldElIds.Print();
//    std::cout << "\n\n";

    elIds[0] = k;
    elIds[1] = k;

    for (int i=1; i<m_numLevels; i++)
    {
        if (elIds[i] == oldElIds[i]) {
            break;
        }
        elIds[i+1]
                = hierMeshTrans[m_numLevels-1-i]->getParentId(elIds[i]);
    }
}


std::tuple<Vector, Vector, Vector>
sparseHeat::SpatialErrorOfSolutionInNaturalNorm
:: eval (double t,
         Array<double> &temporalBasisVals,
         mfem::Array<double> &temporalBasisGradientVals,
         Array<std::shared_ptr<GridFunction>> &spatialTemperatures,
         Array<std::shared_ptr<GridFunction>> &spatialHeatFluxes) const
{
    assert(temporalBasisVals.Size() == m_numLevels+1);
    assert(spatialTemperatures.Size() == temporalBasisVals.Size());
    assert(spatialHeatFluxes.Size() == temporalBasisVals.Size());

    setCurrentTimeForCoefficients(t);

    auto meshes = m_spatialMeshHierarchy->getMeshes();
    auto finestSpatialMesh = meshes[m_numLevels-1];

    auto finestFesForTemperature = spatialTemperatures[0]->FESpace();
//    auto finestFesForHeatFlux = spatialHeatFluxes[0]->FESpace();

    int dim = finestSpatialMesh->Dimension();
    Vector coords(finestSpatialMesh->Dimension());
    IntegrationPoint ipCoarse;

    Array<int> parentElIds(m_numLevels+1);
    initializeSpatialParentElementIds(0, parentElIds);

    Vector temperatureRelErrorNormH1(2);
    temperatureRelErrorNormH1 = 0.;

    Vector heatFluxRelErrorNormL2(2);
    heatFluxRelErrorNormL2 = 0.;

    Vector solDivergenceRelErrorNormL2(2);
    solDivergenceRelErrorNormL2 = 0.;

    for (int k=0; k<finestSpatialMesh->GetNE(); k++)
    {
        auto elTransFine = finestSpatialMesh->GetElementTransformation(k);
        auto feFineForTemperature = finestFesForTemperature->GetFE(k);

        getSpatialParentElementIds(k, parentElIds);

        // integration rules
        int order = 2*feFineForTemperature->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(feFineForTemperature->GetGeomType(),
                                order);

        for (int i=0; i<ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ipFine = ir->IntPoint(i);
            elTransFine->SetIntPoint(&ipFine);
            elTransFine->Transform(ipFine, coords);

            double weight = ipFine.weight * elTransFine->Weight();

            // exact and discrete temperature values
            double uRef
                    = m_exactTemperatureCoeff->Eval(*elTransFine, ipFine);

            double u = temporalBasisVals[0]
                    *spatialTemperatures[0]
                    ->GetValue(k, ipFine);
            u += temporalBasisVals[1]
                    *spatialTemperatures[1]
                    ->GetValue(k, ipFine);
            for (int m=2; m<m_numLevels+1; m++)
            {
                auto elTransCoarse
                        = meshes[m_numLevels-m]
                        ->GetElementTransformation(parentElIds[m]);
                elTransCoarse->TransformBack(coords, ipCoarse);

                u += temporalBasisVals[m]
                        *spatialTemperatures[m]
                        ->GetValue(parentElIds[m], ipCoarse);
            }
//            std::cout << k << "\t" << i << "\t"
//                      << uRef << "\t" << u << std::endl;

            // exact and discrete temperature gradient values
            Vector bufErrorGradu(dim);
            bufErrorGradu = 0.;
            m_exactSpatialGradientOfTemperatureCoeff
                    ->Eval(bufErrorGradu, *elTransFine, ipFine);

            Vector gradu(dim); gradu = 0.;
            Vector graduBuf(dim);
            spatialTemperatures[0]->GetGradient(*elTransFine, graduBuf);
            gradu.Add(temporalBasisVals[0], graduBuf);
            spatialTemperatures[1]->GetGradient(*elTransFine, graduBuf);
            gradu.Add(temporalBasisVals[1], graduBuf);
            for (int m=2; m<m_numLevels+1; m++)
            {
                auto elTransCoarse
                        = meshes[m_numLevels-m]
                        ->GetElementTransformation(parentElIds[m]);
                elTransCoarse->TransformBack(coords, ipCoarse);
                elTransCoarse->SetIntPoint(&ipCoarse);

                spatialTemperatures[m]
                        ->GetGradient(*elTransCoarse, graduBuf);
                gradu.Add(temporalBasisVals[m], graduBuf);
            }

            // discrete temperature temporal gradient values
            double dudt = temporalBasisGradientVals[0]
                    *spatialTemperatures[0]
                    ->GetValue(k, ipFine);
            dudt += temporalBasisGradientVals[1]
                    *spatialTemperatures[1]
                    ->GetValue(k, ipFine);
            for (int m=2; m<m_numLevels+1; m++)
            {
                auto elTransCoarse
                        = meshes[m_numLevels-m]
                        ->GetElementTransformation(parentElIds[m]);
                elTransCoarse->TransformBack(coords, ipCoarse);

                dudt += temporalBasisGradientVals[m]
                        *spatialTemperatures[m]
                        ->GetValue(parentElIds[m], ipCoarse);
            }

            // exact and discrete heat flux values
            Vector bufErrorHeatFlux(dim);
            bufErrorHeatFlux = 0.;
            m_exactHeatFluxCoeff
                    ->Eval(bufErrorHeatFlux, *elTransFine, ipFine);

            Vector heatFlux(dim); heatFlux = 0.;
            Vector heatFluxBuf(dim);
            spatialHeatFluxes[0]->GetVectorValue(k, ipFine, heatFluxBuf);
            heatFlux.Add(temporalBasisVals[0], heatFluxBuf);
            spatialHeatFluxes[1]->GetVectorValue(k, ipFine, heatFluxBuf);
            heatFlux.Add(temporalBasisVals[1], heatFluxBuf);
            for (int m=2; m<m_numLevels+1; m++)
            {
                auto elTransCoarse
                        = meshes[m_numLevels-m]
                        ->GetElementTransformation(parentElIds[m]);
                elTransCoarse->TransformBack(coords, ipCoarse);

                spatialHeatFluxes[m]->GetVectorValue(parentElIds[m], ipCoarse, heatFluxBuf);
                heatFlux.Add(temporalBasisVals[m], heatFluxBuf);
            }

            // discrete heat flux divergence values
            double divxQ = temporalBasisVals[0]
                    *spatialHeatFluxes[0]
                    ->GetDivergence(*elTransFine);
            divxQ += temporalBasisVals[1]
                    *spatialHeatFluxes[1]
                    ->GetDivergence(*elTransFine);
            for (int m=2; m<m_numLevels+1; m++)
            {
                auto elTransCoarse
                        = meshes[m_numLevels-m]
                        ->GetElementTransformation(parentElIds[m]);
                elTransCoarse->TransformBack(coords, ipCoarse);
                elTransCoarse->SetIntPoint(&ipCoarse);

                divxQ += temporalBasisVals[m]
                        *spatialHeatFluxes[m]
                        ->GetDivergence(*elTransCoarse);
            }

            // temperature error
            temperatureRelErrorNormH1(1) += weight
                    * (uRef*uRef + (bufErrorGradu*bufErrorGradu));
            bufErrorGradu.Add(-1, gradu);
            temperatureRelErrorNormH1(0) += weight
                    * ((uRef - u)*(uRef - u)
                       + (bufErrorGradu * bufErrorGradu));

            // heat flux error
            heatFluxRelErrorNormL2(1) += weight
                    * (bufErrorHeatFlux*bufErrorHeatFlux);
            bufErrorHeatFlux.Add(-1, heatFlux);
            heatFluxRelErrorNormL2(0) += weight
                    * (bufErrorHeatFlux*bufErrorHeatFlux);

            // divergence error
            double tmp = m_sourceCoeff->Eval(*elTransFine, ipFine);
            solDivergenceRelErrorNormL2(0) += weight
                    * ((tmp - (dudt - divxQ))*(tmp - (dudt - divxQ)));
            solDivergenceRelErrorNormL2(1) += weight * tmp * tmp;
        }
    }

    return {temperatureRelErrorNormH1,
                heatFluxRelErrorNormL2, solDivergenceRelErrorNormL2};
}

// End of file
