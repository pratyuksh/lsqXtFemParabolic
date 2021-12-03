#ifndef SPARSE_HEAT_SOLUTION_HANDLER_HPP
#define SPARSE_HEAT_SOLUTION_HANDLER_HPP

#include "mfem.hpp"

#include "../mymfem/nested_hierarchy.hpp"


namespace sparseHeat
{

class SolutionHandler
{
public:
    SolutionHandler
    (int minLevel, int maxLevel,
     std::shared_ptr<mymfem::NestedFEHierarchy>&
     spatialNestedFEHierarchyTemperature,
     std::shared_ptr<mymfem::NestedFEHierarchy>&
     spatialNestedFEHierarchyHeatFlux);

    mfem::Vector getTemperatureDataAtInitialTime() const;

    mfem::Vector getTemperatureDataAtEndTime() const;

    mfem::Vector getHeatFluxDataAtInitialTime() const;

    mfem::Vector getHeatFluxDataAtEndTime() const;

    const mfem::Vector& getTemperatureData() const {
        return m_data->GetBlock(0);
    }

    mfem::Vector& getTemperatureData() {
        return m_data->GetBlock(0);
    }

    int getTemperatureDataSize() const {
        return (m_data->GetBlock(0)).Size();
    }

    const mfem::Vector& getHeatFluxData() const {
        return m_data->GetBlock(1);
    }

    mfem::Vector& getHeatFluxData() {
        return m_data->GetBlock(1);
    }

    int getHeatFluxDataSize() const {
        return (m_data->GetBlock(1)).Size();
    }

    std::shared_ptr<mfem::BlockVector> getData() const {
        return m_data;
    }

    mfem::Array<int> getDataSize() const {
        mfem::Array<int> blockSizes(2);
        blockSizes[0] = (m_data->GetBlock(0)).Size();
        blockSizes[1] = (m_data->GetBlock(1)).Size();
        return blockSizes;
    }

private:
    std::shared_ptr<mfem::BlockVector> m_data;

    int m_numLevels;
    mfem::Array<int> m_temporalHierarchicalFESizes;
    mfem::Array<int> m_spatialFESizesTemperature;
    mfem::Array<int> m_spatialFESizesHeatFlux;
};

}

#endif

