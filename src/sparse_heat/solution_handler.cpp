#include "solution_handler.hpp"

#include "utilities.hpp"


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
}

Vector sparseHeat::SolutionHandler
:: getTemperatureDataAtEndTime()
{
     auto temperatureData = getTemperatureData();

     int sizeTemporalCoarsest = m_temporalHierarchicalFESizes[0];
     int sizeSpatialFinest = m_spatialFESizesTemperature[m_numLevels-1];

     Vector temperatureAtEndTime(temperatureData.GetData()
                                 + (sizeTemporalCoarsest-1)*sizeSpatialFinest,
                                 sizeSpatialFinest);

     return temperatureAtEndTime;
}

Vector sparseHeat::SolutionHandler
:: getHeatFluxDataAtEndTime()
{
     auto heatFluxData = getHeatFluxData();

     int sizeTemporalCoarsest = m_temporalHierarchicalFESizes[0];
     int sizeSpatialFinest = m_spatialFESizesHeatFlux[m_numLevels-1];

     Vector heatFluxAtEndTime(heatFluxData.GetData()
                                 + (sizeTemporalCoarsest-1)*sizeSpatialFinest,
                                 sizeSpatialFinest);

     return heatFluxAtEndTime;
}
