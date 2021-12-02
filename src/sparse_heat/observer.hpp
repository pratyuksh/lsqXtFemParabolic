#ifndef SPARSE_HEAT_OBSERVER_HPP
#define SPARSE_HEAT_OBSERVER_HPP

#include "mfem.hpp"

#include "../mymfem/base_observer.hpp"
#include "../mymfem/nested_hierarchy.hpp"
#include "../heat/coefficients.hpp"
#include "discretisation.hpp"
#include "solution_handler.hpp"
#include "spatial_error_evaluator.hpp"


namespace sparseHeat {

/**
 * @brief Observer for the Sparse solver for the Heat equation
 *
 * Provides functions to visualize solutions, compute error
 * and write output to files.
 */
class Observer : public mymfem::BaseObserver
{
public:
    Observer () : BaseObserver () {}

    Observer (const nlohmann::json& config,
              int numLevels,
              int minTemporalLevel);

    //! Sets the test case and discretisation
    void set (std::shared_ptr<heat::TestCases>&,
              std::shared_ptr<LsqSparseXtFem>&);

    void visualizeSolutionAtEndTime(const SolutionHandler&);

    mfem::Vector evalError(const SolutionHandler&);

    /**
     * @brief Evaluates error in the natural norm
     * @param solution handler
     * @return discrete solution error in natural norm
     */
    mfem::Vector evalErrorInNaturalNorm (const SolutionHandler&) const;

private:
    void evalTemporalBasisVals
    (const mfem::Array<double> &temporalMeshSizes,
     const mfem::Array<int> &temporalElIds,
     const mfem::Array<double> &temporalElLeftEdgeCoords,
     double t, mfem::Array<double> &temporalBasisVals) const;

    void evalTemporalBasisGradientVals
    (const mfem::Array<double> &temporalMeshSizes,
     const mfem::Array<int> &temporalElIds,
     const mfem::Array<double> &temporalElLeftEdgeCoords,
     double t, mfem::Array<double> &temporalBasisGradientVals) const;

    void collectSpatialGridFunctionsForTemperature
    (const SolutionHandler&,
     const mfem::Array<int>&,
     const mfem::Array<int>&,
     mfem::Array<std::shared_ptr<mfem::GridFunction>>&) const;

    void collectSpatialGridFunctionsForHeatFlux
    (const SolutionHandler&,
     const mfem::Array<int>&,
     const mfem::Array<int>&,
     mfem::Array<std::shared_ptr<mfem::GridFunction>>&) const;

private:
    int m_numLevels;
    int m_minTemporalLevel, m_maxTemporalLevel;

    double m_endTime;

    bool m_boolEvalError;
    std::string m_errorType;

    std::shared_ptr<heat::TestCases> m_testCase;
    std::shared_ptr<LsqSparseXtFem> m_disc;

    std::shared_ptr<mymfem::NestedMeshHierarchy> m_spatialMeshHierarchy;

    mymfem::NestedFESpaces m_spatialNestedFESpacesForTemperature;
    mymfem::NestedFESpaces m_spatialNestedFESpacesForHeatFlux;

    std::shared_ptr<mfem::FiniteElementSpace> m_spatialFinestFESpaceForTemperature;
    std::shared_ptr<mfem::FiniteElementSpace> m_spatialFinestFESpaceForHeatFlux;

    mutable std::shared_ptr<mfem::MatrixCoefficient> m_materialCoeff;

    mutable std::shared_ptr<mfem::Coefficient> m_exactTemperatureCoeff;
    mutable std::shared_ptr<mfem::VectorCoefficient> m_exactHeatFluxCoeff;

    mutable std::shared_ptr<mfem::VectorCoefficient>
    m_exactSpatialGradientOfTemperatureCoeff;
    mutable std::shared_ptr<mfem::Coefficient>
    m_exactTemporalGradientOfTemperatureCoeff;

    mutable std::shared_ptr<mfem::Coefficient> m_sourceCoeff;

    mutable std::unique_ptr<SpatialErrorOfSolution>
    m_spatialErrorOfSolutionInNaturalNorm;
};

}

#endif // SPARSE_HEAT_OBSERVER_HPP
