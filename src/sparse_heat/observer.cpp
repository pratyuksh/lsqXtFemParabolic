#include "observer.hpp"

#include "spatial_error_evaluator.hpp"
#include "utilities.hpp"

using namespace mfem;


sparseHeat::Observer
:: Observer (const nlohmann::json& config,
             int numLevels,
             int minTemporalLevel)
    : BaseObserver (config),
      m_numLevels(numLevels),
      m_minTemporalLevel(minTemporalLevel)
{
    m_maxTemporalLevel = m_minTemporalLevel + m_numLevels - 1;

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "end_time",
                                        m_endTime, 1.);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "eval_error",
                                        m_boolEvalError, false);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "error_type",
                                        m_errorType, "natural");
}

void sparseHeat::Observer
:: set (std::shared_ptr<heat::TestCases> &testCase,
        std::shared_ptr<sparseHeat::LsqSparseXtFem> &disc)
{
    m_testCase = testCase;
    m_disc = disc;

    m_spatialMeshHierarchy
            = disc->getSpatialNestedFEHierarchyForTemperature()
            ->getNestedMeshHierarchy();
    assert(m_spatialMeshHierarchy->getNumMeshes() == m_numLevels);

    m_spatialNestedFESpacesForTemperature
            = disc->getSpatialNestedFEHierarchyForTemperature()
            ->getFESpaces();
    m_spatialFinestFESpaceForTemperature
            = m_spatialNestedFESpacesForTemperature[m_numLevels-1];

    m_spatialNestedFESpacesForHeatFlux
            = disc->getSpatialNestedFEHierarchyForHeatFlux()
            ->getFESpaces();
    m_spatialFinestFESpaceForHeatFlux
            = m_spatialNestedFESpacesForHeatFlux[m_numLevels-1];

    if (m_boolEvalError)
    {
        m_materialCoeff = std::make_shared<heat::MediumTensorCoeff>
                (m_testCase);

        m_exactTemperatureCoeff
                = std::make_shared<heat::ExactTemperatureCoeff>(m_testCase);
        m_exactHeatFluxCoeff
                = std::make_shared<heat::ExactHeatFluxCoeff>(m_testCase);

        m_exactSpatialGradientOfTemperatureCoeff
                = std::make_shared<heat::ExactTemperatureSpatialGradCoeff>
                (m_testCase);
        m_exactTemporalGradientOfTemperatureCoeff
                = std::make_shared<heat::ExactTemperatureTemporalGradCoeff>
                (m_testCase);

        m_sourceCoeff = std::make_shared<heat::SourceCoeff>(m_testCase);

        m_spatialErrorOfSolutionInNaturalNorm
                    = std::make_unique<SpatialErrorOfSolutionInNaturalNorm>
                    (m_spatialMeshHierarchy);
        m_spatialErrorOfSolutionInNaturalNorm
                ->setCoefficients(m_materialCoeff,
                                  m_exactTemperatureCoeff,
                                  m_exactHeatFluxCoeff,
                                  m_exactSpatialGradientOfTemperatureCoeff,
                                  m_exactTemporalGradientOfTemperatureCoeff,
                                  m_sourceCoeff);

        m_spatialErrorOfSolutionInLeastSquaresNorm
                    = std::make_unique<SpatialErrorOfSolutionInLeastSquaresNorm>
                    (m_spatialMeshHierarchy);
        m_spatialErrorOfSolutionInLeastSquaresNorm
                ->setCoefficients(m_materialCoeff,
                                  m_exactTemperatureCoeff,
                                  m_exactHeatFluxCoeff,
                                  m_exactSpatialGradientOfTemperatureCoeff,
                                  m_exactTemporalGradientOfTemperatureCoeff,
                                  m_sourceCoeff);
    }
}

void sparseHeat::Observer
:: visualizeSolutionAtEndTime
(const sparseHeat::SolutionHandler& solutionHandler)
{
    auto spatialFESpacesForTemperature
            = m_disc->getSpatialNestedFEHierarchyForTemperature()
            ->getFESpaces();
    auto temperatureDataAtEndTime
            = solutionHandler.getTemperatureDataAtEndTime();

    auto temperatureAtEndTime
            = std::make_unique<GridFunction>
            (spatialFESpacesForTemperature[m_numLevels-1].get(),
            temperatureDataAtEndTime.GetData());
    visualize(*temperatureAtEndTime);

    auto spatialFESpacesForHeatFlux
            = m_disc->getSpatialNestedFEHierarchyForHeatFlux()
            ->getFESpaces();
    auto heatFluxDataAtEndTime
            = solutionHandler.getHeatFluxDataAtEndTime();

    auto heatFluxAtEndTime
            = std::make_unique<GridFunction>
            (spatialFESpacesForHeatFlux[m_numLevels-1].get(),
            heatFluxDataAtEndTime.GetData());
    visualize(*heatFluxAtEndTime);
}

Vector sparseHeat::Observer
:: evalError(const sparseHeat::SolutionHandler& solutionHandler)
{
    Vector solutionError;

    if (m_errorType == "natural") {
        solutionError = evalErrorInNaturalNorm(solutionHandler);
    }
    else if (m_errorType == "lsq") {
        solutionError = evalErrorInLeastSquaresNorm(solutionHandler);
    }

    return solutionError;
}

mfem::Vector sparseHeat::Observer
:: evalErrorInNaturalNorm(const SolutionHandler & solutionHandler) const
{
    Vector solutionError(3);

    Vector temperatureRelErrorNormL2H1(2);
    Vector heatFluxRelErrorNormL2L2(2);
    Vector solDivergenceRelErrorNormL2L2(2);

    // temporal integration rule
    int order = 3;
    const IntegrationRule *ir
            = &IntRules.Get(Geometry::SEGMENT, order);

    Array<double> temporalMeshSizes
            = evalTemporalMeshSizes(m_endTime,
                                    m_minTemporalLevel,
                                    m_numLevels);
    Array<int> temporalFesDims
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);

    // auxiliary variables
    Array<int> temporalElIds(m_numLevels);
    Array<double> temporalElLeftEdgeCoords(m_numLevels);

    Array<double> temporalBasisVals(m_numLevels+1);
    Array<double> temporalBasisGradientVals(m_numLevels+1);
    Array<std::shared_ptr<GridFunction>> spatialTemperatures(m_numLevels+1);
    Array<std::shared_ptr<GridFunction>> spatialHeatFluxes(m_numLevels+1);

    int numElsTemporalFinestMesh
            = static_cast<int>(std::pow(2, m_maxTemporalLevel));
    double finestTemporalMeshSize = temporalMeshSizes[m_numLevels-1];
//    std::cout << "\nNumber of temporal elements on finest level: "
//              << numElsTemporalFinestMesh << std::endl;

    temperatureRelErrorNormL2H1 = 0.;
    heatFluxRelErrorNormL2L2 = 0.;
    solDivergenceRelErrorNormL2L2 = 0.;
    for (int n=0; n<numElsTemporalFinestMesh; n++)
    {
        // find current temporal element id for all temporal mesh levels
        evalTemporalElIds(n, temporalElIds);

        // evaluate the coordinate of left edge
        // of the current temporal element for all temporal mesh levels
        evalTemporalElLeftEdgeCoords(temporalMeshSizes, temporalElIds,
                                     temporalElLeftEdgeCoords);

        collectSpatialGridFunctionsForTemperature
                (solutionHandler, temporalFesDims, temporalElIds,
                 spatialTemperatures);
        collectSpatialGridFunctionsForHeatFlux
                (solutionHandler, temporalFesDims, temporalElIds,
                 spatialHeatFluxes);

        double tFinestLeft = temporalElLeftEdgeCoords[m_numLevels-1];
        double tFinestRight = tFinestLeft + temporalMeshSizes[m_numLevels-1];

//        std::cout << "\n" << tFinestLeft << "\t"
//                  << tFinestRight << "\n\n";

        Vector localTemperatureRelErrorNormH1;
        Vector localHeatFluxRelErrorNormL2;
        Vector localSolDivergenceRelErrorNormL2;
        for (int i=0; i<ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            double t = affineTransform(tFinestLeft, tFinestRight, ip.x);

            evalTemporalBasisVals(temporalMeshSizes,
                                  temporalElIds,
                                  temporalElLeftEdgeCoords,
                                  t, temporalBasisVals);
            evalTemporalBasisGradientVals(temporalMeshSizes,
                                          temporalElIds,
                                          temporalElLeftEdgeCoords,
                                          t, temporalBasisGradientVals);

            std::tie (localTemperatureRelErrorNormH1,
                      localHeatFluxRelErrorNormL2,
                      localSolDivergenceRelErrorNormL2)
                    = m_spatialErrorOfSolutionInNaturalNorm
                    ->eval(t, temporalBasisVals, temporalBasisGradientVals,
                           spatialTemperatures, spatialHeatFluxes);

            temperatureRelErrorNormL2H1
                    .Add(ip.weight*finestTemporalMeshSize,
                         localTemperatureRelErrorNormH1);
            heatFluxRelErrorNormL2L2
                    .Add(ip.weight*finestTemporalMeshSize,
                         localHeatFluxRelErrorNormL2);
            solDivergenceRelErrorNormL2L2
                    .Add(ip.weight*finestTemporalMeshSize,
                         localSolDivergenceRelErrorNormL2);
        }
    }
    solutionError(0) = std::sqrt(temperatureRelErrorNormL2H1(0))
            / std::sqrt(temperatureRelErrorNormL2H1(1));
    solutionError(1) = std::sqrt(heatFluxRelErrorNormL2L2(0))
            / std::sqrt(heatFluxRelErrorNormL2L2(1));
    solutionError(2) = std::sqrt(solDivergenceRelErrorNormL2L2(0))
            / std::sqrt(solDivergenceRelErrorNormL2L2(1));

    return solutionError;
}

mfem::Vector sparseHeat::Observer
:: evalErrorInLeastSquaresNorm(const SolutionHandler & solutionHandler) const
{
    Vector solutionError(3);

    Vector pdeErrorNormL2L2(1);
    Vector heatFluxErrorNormL2L2(1);
    Vector initialTemperatureErrorNormL2(1);

    // temporal integration rule
    int order = 3;
    const IntegrationRule *ir
            = &IntRules.Get(Geometry::SEGMENT, order);

    Array<double> temporalMeshSizes
            = evalTemporalMeshSizes(m_endTime,
                                    m_minTemporalLevel,
                                    m_numLevels);
    Array<int> temporalFesDims
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);

    // auxiliary variables
    Array<int> temporalElIds(m_numLevels);
    Array<double> temporalElLeftEdgeCoords(m_numLevels);

    Array<double> temporalBasisVals(m_numLevels+1);
    Array<double> temporalBasisGradientVals(m_numLevels+1);
    Array<std::shared_ptr<GridFunction>> spatialTemperatures(m_numLevels+1);
    Array<std::shared_ptr<GridFunction>> spatialHeatFluxes(m_numLevels+1);

    int numElsTemporalFinestMesh
            = static_cast<int>(std::pow(2, m_maxTemporalLevel));
    double finestTemporalMeshSize = temporalMeshSizes[m_numLevels-1];

    // pde and heat flux error
    pdeErrorNormL2L2 = 0.;
    heatFluxErrorNormL2L2 = 0.;
    for (int n=0; n<numElsTemporalFinestMesh; n++)
    {
        // find current temporal element id for all temporal mesh levels
        evalTemporalElIds(n, temporalElIds);

        // evaluate the coordinate of left edge
        // of the current temporal element for all temporal mesh levels
        evalTemporalElLeftEdgeCoords(temporalMeshSizes, temporalElIds,
                                     temporalElLeftEdgeCoords);

        collectSpatialGridFunctionsForTemperature
                (solutionHandler, temporalFesDims, temporalElIds,
                 spatialTemperatures);
        collectSpatialGridFunctionsForHeatFlux
                (solutionHandler, temporalFesDims, temporalElIds,
                 spatialHeatFluxes);

        double tFinestLeft = temporalElLeftEdgeCoords[m_numLevels-1];
        double tFinestRight = tFinestLeft + temporalMeshSizes[m_numLevels-1];

        Vector localPdeErrorNormL2;
        Vector localHeatFluxErrorNormL2;
        for (int i=0; i<ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            double t = affineTransform(tFinestLeft, tFinestRight, ip.x);

            evalTemporalBasisVals(temporalMeshSizes,
                                  temporalElIds,
                                  temporalElLeftEdgeCoords,
                                  t, temporalBasisVals);
            evalTemporalBasisGradientVals(temporalMeshSizes,
                                          temporalElIds,
                                          temporalElLeftEdgeCoords,
                                          t, temporalBasisGradientVals);

            std::tie (localPdeErrorNormL2,
                      localHeatFluxErrorNormL2,
                      std::ignore)
                    = m_spatialErrorOfSolutionInLeastSquaresNorm
                    ->eval(t, temporalBasisVals, temporalBasisGradientVals,
                           spatialTemperatures, spatialHeatFluxes);

            pdeErrorNormL2L2
                    .Add(ip.weight*finestTemporalMeshSize,
                         localPdeErrorNormL2);
            heatFluxErrorNormL2L2
                    .Add(ip.weight*finestTemporalMeshSize,
                         localHeatFluxErrorNormL2);
        }
    }
    solutionError(0) = std::sqrt(pdeErrorNormL2L2(0));
    solutionError(1) = std::sqrt(heatFluxErrorNormL2L2(0));

    // error in initial conditions
    auto initialTemperatureData
            = solutionHandler.getTemperatureDataAtInitialTime();
    auto initialTemperature
            = std::make_unique<GridFunction>
            (m_spatialFinestFESpaceForTemperature.get(),
             initialTemperatureData.GetData());
    m_exactTemperatureCoeff->SetTime(0.);
    solutionError(2)
            = initialTemperature->ComputeL2Error(*m_exactTemperatureCoeff);

    return solutionError;
}

void sparseHeat::Observer
:: evalTemporalBasisVals
(const mfem::Array<double> &temporalMeshSizes,
 const mfem::Array<int> &temporalElIds,
 const mfem::Array<double> &temporalElLeftEdgeCoords,
 double t,
 mfem::Array<double> &temporalBasisVals) const
{
    // coarsest temporal level has standard FE basis functions
    {
        temporalBasisVals[0]
                = evalRightHalfOfHatBasis
                (t, temporalElLeftEdgeCoords[0],
                 temporalMeshSizes[0]);

        temporalBasisVals[1]
                = evalLeftHalfOfHatBasis
                (t, temporalElLeftEdgeCoords[0],
                 temporalMeshSizes[0]);
    }
//    std::cout << temporalBasisVals[0] << "\t"
//              << temporalBasisVals[1] << std::endl;
    // remaining temporal levels have hierarchical basis functions
    for (int m=1; m<m_numLevels; m++)
    {
        if (temporalElIds[m]%2 == 0) {
            temporalBasisVals[m+1]
                    = evalLeftHalfOfHatBasis
                    (t, temporalElLeftEdgeCoords[m],
                     temporalMeshSizes[m]);
        }
        else {
            temporalBasisVals[m+1]
                    = evalRightHalfOfHatBasis
                    (t, temporalElLeftEdgeCoords[m],
                     temporalMeshSizes[m]);
        }
//        std::cout << temporalBasisVals[m+1] << std::endl;
    }
}

void sparseHeat::Observer
:: evalTemporalBasisGradientVals
(const mfem::Array<double> &temporalMeshSizes,
 const mfem::Array<int> &temporalElIds,
 const mfem::Array<double> &temporalElLeftEdgeCoords,
 double t,
 mfem::Array<double> &temporalBasisGradientVals) const
{
    // coarsest temporal level has standard FE basis functions
    {
        temporalBasisGradientVals[0] = -1./temporalMeshSizes[0];
        temporalBasisGradientVals[1] = +1./temporalMeshSizes[0];
    }
    // remaining temporal levels have hierarchical basis functions
    for (int m=1; m<m_numLevels; m++)
    {
        if (temporalElIds[m]%2 == 0) {
            temporalBasisGradientVals[m+1] = +1./temporalMeshSizes[m];
        }
        else {
            temporalBasisGradientVals[m+1] = -1./temporalMeshSizes[m];
        }
    }
}

void sparseHeat::Observer
:: collectSpatialGridFunctionsForTemperature
(const SolutionHandler& solutionHandler,
 const Array<int>& temporalFesDims,
 const Array<int>& temporalElIds,
 Array<std::shared_ptr<GridFunction>>& spatialTemperature) const
{
    auto solutionData = solutionHandler.getData();

    int shift = 0;

    // coarsest temporal level has standard FE hat basis functions
    // so there are two grid functions at the finest spatial mesh level
    {
        const int m = 0;

        int i = temporalElIds[m];
        int spatialId = m_numLevels - m - 1;

        auto spatialFes
                = m_spatialNestedFESpacesForTemperature[spatialId];
        int spatialFesDim
                = spatialFes->GetTrueVSize();

        spatialTemperature[m]
                = std::make_shared<GridFunction>
                (spatialFes.get(),
                 solutionData->GetData() + i*spatialFesDim + shift);

        spatialTemperature[m+1]
                = std::make_shared<GridFunction>
                (spatialFes.get(),
                 solutionData->GetData() + (i+1)*spatialFesDim + shift);

        shift += temporalFesDims[m] * spatialFesDim;
    }
    // remaining levels
    for (int m=1; m<m_numLevels; m++)
    {
        int i = temporalElIds[m]/2;
        int spatialId = m_numLevels - m - 1;

        auto spatialFes
                = m_spatialNestedFESpacesForTemperature[spatialId];
        int spatialFesDim
                = spatialFes->GetTrueVSize();

        spatialTemperature[m+1]
                = std::make_shared<GridFunction>
                (spatialFes.get(),
                 solutionData->GetData() + i*spatialFesDim + shift);

        shift += temporalFesDims[m] * spatialFesDim;
    }
}

void sparseHeat::Observer
:: collectSpatialGridFunctionsForHeatFlux
(const SolutionHandler& solutionHandler,
 const Array<int>& temporalFesDims,
 const Array<int>& temporalElIds,
 Array<std::shared_ptr<GridFunction>>& spatialHeatFlux) const
{
    auto solutionData = solutionHandler.getData();

    int shift = solutionHandler.getTemperatureDataSize();

    // coarsest temporal level has standard FE hat basis functions
    // so there are two grid functions at the finest spatial mesh level
    {
        const int m = 0;

        int i = temporalElIds[m];
        int spatialId = m_numLevels - m - 1;

        auto spatialFes
                = m_spatialNestedFESpacesForHeatFlux[spatialId];
        int spatialFesDim
                = spatialFes->GetTrueVSize();

        spatialHeatFlux[m]
                = std::make_shared<GridFunction>
                (spatialFes.get(),
                 solutionData->GetData() + i*spatialFesDim + shift);

        spatialHeatFlux[m+1]
                = std::make_shared<GridFunction>
                (spatialFes.get(),
                 solutionData->GetData() + (i+1)*spatialFesDim + shift);

        shift += temporalFesDims[m] * spatialFesDim;
    }
    // remaining levels
    for (int m=1; m<m_numLevels; m++)
    {
        int i = temporalElIds[m]/2;
        int spatialId = m_numLevels - m - 1;

        auto spatialFes
                = m_spatialNestedFESpacesForHeatFlux[spatialId];
        int spatialFesDim
                = spatialFes->GetTrueVSize();

        spatialHeatFlux[m+1]
                = std::make_shared<GridFunction>
                (spatialFes.get(),
                 solutionData->GetData() + i*spatialFesDim + shift);

        shift += temporalFesDims[m] * spatialFesDim;
    }
}
