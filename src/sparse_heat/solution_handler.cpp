#include "solution_handler.hpp"

#include "utilities.hpp"

using namespace mfem;


sparseHeat::SolutionHandler
:: SolutionHandler
(int minLevel, int maxLevel,
 std::shared_ptr<mymfem::NestedFEHierarchy> &spatialNestedFEHierarchyTemperature,
 std::shared_ptr<mymfem::NestedFEHierarchy> &spatialNestedFEHierarchyHeatFlux)
{
    m_numLevels = (maxLevel - minLevel) + 1;
    m_temporalHierarchicalFESizes
            = evalTemporalBlockSizes(minLevel, maxLevel);
    m_spatialFESizesTemperature
            = spatialNestedFEHierarchyTemperature->getNumDims();
    m_spatialFESizesHeatFlux
            = spatialNestedFEHierarchyHeatFlux->getNumDims();

    // size of space-time solution vector
    auto tmpBuf1 = evalSpaceTimeBlockSizes(m_spatialFESizesTemperature,
                                           m_temporalHierarchicalFESizes);
    int temperatureDataSize = tmpBuf1.Sum();

    auto tmpBuf2 = evalSpaceTimeBlockSizes(m_spatialFESizesHeatFlux,
                                           m_temporalHierarchicalFESizes);
    int heatFluxDataSize = tmpBuf2.Sum();

    Array<int> blockOffsets(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = temperatureDataSize;
    blockOffsets[2] = blockOffsets[1] + heatFluxDataSize;

    m_data = std::make_shared<BlockVector>(blockOffsets);
    (*m_data) = 0.;
}

Vector sparseHeat::SolutionHandler
:: getTemperatureDataAtInitialTime() const
{
    auto temperatureData = getTemperatureData();

    int sizeSpatialFinest = m_spatialFESizesTemperature[m_numLevels-1];
    Vector temperatureAtInitialTime(sizeSpatialFinest);
    temperatureAtInitialTime = temperatureData.GetData();

    return temperatureAtInitialTime;
}

Vector sparseHeat::SolutionHandler
:: getTemperatureDataAtEndTime() const
{
    auto temperatureData = getTemperatureData();

    int sizeTemporalCoarsest = m_temporalHierarchicalFESizes[0];
    int sizeSpatialFinest = m_spatialFESizesTemperature[m_numLevels-1];

    int offsets = (sizeTemporalCoarsest-1)*sizeSpatialFinest;
    Vector temperatureAtEndTime(temperatureData.GetData() + offsets,
                                sizeSpatialFinest);

    return temperatureAtEndTime;
}

Vector sparseHeat::SolutionHandler
:: getHeatFluxDataAtInitialTime() const
{
    auto heatFluxData = getHeatFluxData();

    int sizeSpatialFinest = m_spatialFESizesHeatFlux[m_numLevels-1];
    Vector heatFluxAtInitialTime(sizeSpatialFinest);
    heatFluxAtInitialTime = heatFluxData.GetData();

    return heatFluxAtInitialTime;
}

Vector sparseHeat::SolutionHandler
:: getHeatFluxDataAtEndTime() const
{
    auto heatFluxData = getHeatFluxData();

    int sizeTemporalCoarsest = m_temporalHierarchicalFESizes[0];
    int sizeSpatialFinest = m_spatialFESizesHeatFlux[m_numLevels-1];

    int offsets = (sizeTemporalCoarsest-1)*sizeSpatialFinest;
    Vector heatFluxAtEndTime(heatFluxData.GetData() + offsets,
                             sizeSpatialFinest);

    return heatFluxAtEndTime;
}
