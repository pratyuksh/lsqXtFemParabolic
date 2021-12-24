#ifndef SPARSE_HEAT_SOLVER_HPP
#define SPRASE_HEAT_SOLVER_HPP

#include "mfem.hpp"

#include "../core/config.hpp"
#include "../pardiso/pardiso.hpp"

#include "../heat/test_cases_factory.hpp"
#include "../heat/coefficients.hpp"
#include "solution_handler.hpp"
#include "discretisation.hpp"


namespace sparseHeat{

/**
 * @brief Least-sqares sparse space-time FEM solver for the heat equation
 */
class Solver
{
public:
    Solver (const nlohmann::json& config,
            std::shared_ptr<heat::TestCases>& testCase,
            std::string meshDir,
            const int numLevels,
            const int minSpatialLevel,
            const int minTemporalLevel,
            const bool loadInitMesh=false);

    ~ Solver ();

    void setConfigParams();

    void setMeshHierarchy();

    void setDiscretisation();

    void run();
    std::pair<mfem::Vector, int> runAndMeasurePerformanceMetrics();

    void initialize();

    void assembleSystem();
    void assembleRhs();

    void solve();
    std::pair<double, int> solveAndMeasurePerformanceMetrics ();
    std::pair<double, int> solve (const mfem::Vector&, mfem::Vector&);
    void setPardisoSolver();
    void finalizePardisoSolver();

public:
    double getMeshwidthOfFinestTemporalMesh();

    double getMeshwidthOfFinestSpatialMesh();

    int getNumDofs() {
        return m_solutionHandler->getDataSize().Sum();
    }

public:
    std::shared_ptr<heat::TestCases> getTestCase() const {
        return m_testCase;
    }

    std::shared_ptr<sparseHeat::LsqSparseXtFem> getDiscretisation() const {
        return m_disc;
    }

    std::shared_ptr<sparseHeat::SolutionHandler> getSolutionHandler() const {
        return m_solutionHandler;
    }

private:
    const nlohmann::json& m_config;

    std::shared_ptr<heat::TestCases> m_testCase;

    std::string m_meshDir;
    std::string m_meshElemType;

    int m_numLevels;
    int m_minSpatialLevel;
    int m_minTemporalLevel, m_maxTemporalLevel;
    bool m_loadInitMesh;

    std::string m_discType;
    std::string m_linearSolver;

    double m_endTime;

    std::shared_ptr<mymfem::NestedMeshHierarchy> m_spatialMeshHierarchy;

    std::shared_ptr<sparseHeat::LsqSparseXtFem> m_disc;

    std::shared_ptr<sparseHeat::SolutionHandler> m_solutionHandler;
    std::shared_ptr <mfem::BlockVector> m_rhs;

    mfem::SparseMatrix *m_systemMat = nullptr;
//    mfem::BlockOperator *m_systemOp = nullptr;

#ifdef PARDISO_HPP
    std::unique_ptr<PardisoSolver> m_pardisoSolver;
#endif
};

}

#endif // SPARSE_HEAT_SOLVER_HPP
