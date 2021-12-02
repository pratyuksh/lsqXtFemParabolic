#ifndef SPARSE_HEAT_SPATIAL_ERROR_EVALUATOR_HPP
#define SPARSE_HEAT_SPATIAL_ERROR_EVALUATOR_HPP

#include "mfem.hpp"

#include "../mymfem/nested_hierarchy.hpp"

namespace sparseHeat
{

class SpatialErrorOfSolution
{
public:
    SpatialErrorOfSolution
    (std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialMeshHierarchy) :
    m_spatialMeshHierarchy (spatialMeshHierarchy) {
        m_numLevels = m_spatialMeshHierarchy->getNumMeshes();
    }

    void setCoefficients(std::shared_ptr<mfem::Coefficient> &,
                         std::shared_ptr<mfem::VectorCoefficient> &,
                         std::shared_ptr<mfem::VectorCoefficient> &,
                         std::shared_ptr<mfem::Coefficient> &,
                         std::shared_ptr<mfem::Coefficient> &);

protected:
    void setCurrentTimeForCoefficients(double) const;

    void initializeSpatialParentElementIds
    (int k, mfem::Array<int>& elIds) const;

    void getSpatialParentElementIds
    (int k, mfem::Array<int>& elIds) const;

public:
    virtual std::tuple<mfem::Vector, mfem::Vector, mfem::Vector>
    eval (double, mfem::Array<double>&, mfem::Array<double>&,
          mfem::Array<std::shared_ptr<mfem::GridFunction>>&,
          mfem::Array<std::shared_ptr<mfem::GridFunction>>&) const = 0;

protected:
    int m_numLevels;
    std::shared_ptr<mymfem::NestedMeshHierarchy> m_spatialMeshHierarchy;

    mutable std::shared_ptr<mfem::Coefficient> m_exactTemperatureCoeff;
    mutable std::shared_ptr<mfem::VectorCoefficient> m_exactHeatFluxCoeff;

    mutable std::shared_ptr<mfem::VectorCoefficient>
    m_exactSpatialGradientOfTemperatureCoeff;
    mutable std::shared_ptr<mfem::Coefficient>
    m_exactTemporalGradientOfTemperatureCoeff;

    mutable std::shared_ptr<mfem::Coefficient> m_sourceCoeff;
};


class SpatialErrorOfSolutionInNaturalNorm
        : public SpatialErrorOfSolution
{
public:
    SpatialErrorOfSolutionInNaturalNorm
    (std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialMeshHierarchy) :
    SpatialErrorOfSolution(spatialMeshHierarchy) {}

    std::tuple<mfem::Vector, mfem::Vector, mfem::Vector>
    eval (double, mfem::Array<double>&, mfem::Array<double>&,
          mfem::Array<std::shared_ptr<mfem::GridFunction>>&,
          mfem::Array<std::shared_ptr<mfem::GridFunction>>&)
    const override;
};

}

#endif
