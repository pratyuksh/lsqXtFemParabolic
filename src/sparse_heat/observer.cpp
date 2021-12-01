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

    m_endTime = 1;
    if (config.contains("end_time")) {
        m_endTime = config["end_time"];
    }

    m_boolEvalError = false;
    if (config.contains("eval_error")) {
        m_boolEvalError = config["eval_error"];
    }

    m_errorType = "natural";
    if (config.contains("error_type")) {
        m_errorType = config["error_type"];
    }
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
        m_materialCoeff = std::make_unique<heat::MediumTensorCoeff>
                (m_testCase);

        m_exactTemperatureCoeff
                = std::make_unique<heat::ExactTemperatureCoeff>(m_testCase);
        m_exactHeatFluxCoeff
                = std::make_unique<heat::ExactFluxCoeff>(m_testCase);

        m_exactTemporalGradientOfTemperatureCoeff
                = std::make_unique<heat::ExactTemperatureTimeGradCoeff>
                (m_testCase);

        m_exactTemperature
                = std::make_shared<GridFunction>
                (m_spatialFinestFESpaceForTemperature.get());
        m_exactHeatFlux
                = std::make_shared<GridFunction>
                (m_spatialFinestFESpaceForHeatFlux.get());

        m_spatialErrorOfGradientOfTemperature
                    = std::make_unique<SpatialErrorOfGradientOfTemperature>
                    (m_spatialMeshHierarchy);
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

    evalErrorInNaturalNorm(solutionHandler);

    return solutionError;
}

void sparseHeat::Observer
:: evalErrorInNaturalNorm(const SolutionHandler & solutionHandler) const
{
    double temperatureErrorNormL2H1 = 0;
    double exactTemperatureNormL2H1 = 0;

    // temporal integration rule
    int order = 4;
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

    Array<std::shared_ptr<GridFunction>> spatialTemperature(m_numLevels+1);
    Array<std::shared_ptr<GridFunction>> spatialHeatFlux(m_numLevels+1);

    int numElsTemporalFinestMesh = temporalMeshSizes[m_numLevels-1];

    // buffer variables
    Array<double> bufTemperatureErrorNormL2H1(numElsTemporalFinestMesh);
    bufTemperatureErrorNormL2H1 = 0.;

    Array<double> bufExactTemperatureNormL2H1(numElsTemporalFinestMesh);
    bufExactTemperatureNormL2H1 = 0.;

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
                 spatialTemperature);
        collectSpatialGridFunctionsForHeatFlux
                (solutionHandler, temporalFesDims, temporalElIds,
                 spatialHeatFlux);

        double tFinestLeft = temporalElLeftEdgeCoords[m_numLevels-1];
        double tFinestRight = tFinestLeft + temporalMeshSizes[m_numLevels-1];

        for (int i=0; i<ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            double t = affineTransform(tFinestLeft, tFinestRight, ip.x);

            m_exactTemperatureCoeff->SetTime(t);
            m_exactTemperature
                    ->ProjectCoefficient(*m_exactTemperatureCoeff);

            evalTemporalBasisVals(temporalMeshSizes,
                                  temporalElIds,
                                  temporalElLeftEdgeCoords,
                                  t, temporalBasisVals);

            double tmpVal
                    = m_spatialErrorOfGradientOfTemperature
                    ->eval(m_exactTemperature,
                           temporalBasisVals,
                           spatialTemperature);
        }
    }
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
                (t, temporalElLeftEdgeCoords[1],
                 temporalMeshSizes[1]);
    }
    // remaining temporal levels have hierarchical basis functions
    for (int m=2; m<m_numLevels+1; m++)
    {
        if (temporalElIds[m]%2 == 0) {
            temporalBasisVals[m]
                    = evalLeftHalfOfHatBasis
                    (t, temporalElLeftEdgeCoords[m],
                     temporalMeshSizes[m]);
        }
        else {
            temporalBasisVals[m]
                    = evalRightHalfOfHatBasis
                    (t, temporalElLeftEdgeCoords[m],
                     temporalMeshSizes[m]);
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
{}
