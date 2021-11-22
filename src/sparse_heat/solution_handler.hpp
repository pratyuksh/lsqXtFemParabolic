#ifndef SPARSE_HEAT_SOLUTION_HANDLER_HPP
#define SPARSE_HEAT_SOLUTION_HANDLER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../mymfem/nested_hierarchy.hpp"


namespace sparseHeat
{

class SolutionHandler
{
public:
    SolutionHandler
    (int minLevel, int maxLevel,
     std::shared_ptr<mymfem::NestedFEHierarchy>& spatialNestedFEHierarchyTemperature,
     std::shared_ptr<mymfem::NestedFEHierarchy>& spatialNestedFEHierarchyHeatFlux);

    Vector getTemperatureDataAtEndTime();

    Vector getHeatFluxDataAtEndTime();

    Vector& getTemperatureData() {
        return m_data->GetBlock(0);
    }

    int getTemperatureDataSize() {
        return (m_data->GetBlock(0)).Size();
    }

    Vector& getHeatFluxData() {
        return m_data->GetBlock(1);
    }

    int getHeatFluxDataSize() {
        return (m_data->GetBlock(1)).Size();
    }

    std::shared_ptr<BlockVector> getData() const {
        return m_data;
    }

    Array<int> getDataSize() {
        Array<int> blockSizes(2);
        blockSizes[0] = (m_data->GetBlock(0)).Size();
        blockSizes[1] = (m_data->GetBlock(1)).Size();
        return blockSizes;
    }

private:
    std::shared_ptr<BlockVector> m_data;

    int m_numLevels;
    Array<int> m_temporalHierarchicalFESizes;
    Array<int> m_spatialFESizesTemperature;
    Array<int> m_spatialFESizesHeatFlux;
};

}

#endif

