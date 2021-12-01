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
    m_spatialMeshHierarchy (spatialMeshHierarchy) {}

    virtual double eval
    (std::shared_ptr<mfem::GridFunction>&,
     mfem::Array<double>&,
     mfem::Array<std::shared_ptr<mfem::GridFunction>>&) const = 0;

protected:
    std::shared_ptr<mymfem::NestedMeshHierarchy> m_spatialMeshHierarchy;
};


class SpatialErrorOfGradientOfTemperature : public SpatialErrorOfSolution
{
public:
    SpatialErrorOfGradientOfTemperature
    (std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialMeshHierarchy) :
    SpatialErrorOfSolution(spatialMeshHierarchy) {}

    double eval
    (std::shared_ptr<mfem::GridFunction>&,
     mfem::Array<double>&,
     mfem::Array<std::shared_ptr<mfem::GridFunction>>&) const override;
};

}

#endif
