#ifndef SPARSE_HEAT_SOLVER_HPP
#define SPRASE_HEAT_SOLVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include <Eigen/Core>

#include "../core/config.hpp"
#include "../pardiso/pardiso.hpp"

#include "../heat/test_cases_factory.hpp"
#include "../heat/coefficients.hpp"
#include "solution_handler.hpp"
#include "discretisation.hpp"
#include "observer.hpp"


namespace sparseHeat{

/**
 * @brief Least-sqares sparse space-time FEM solver for the heat equation
 */
class Solver
{
public:
    //! Constructor with config JSON, test case and space-time mesh info
    Solver (const nlohmann::json& config,
            std::shared_ptr<heat::TestCases>& testCase,
            std::string meshDir,
            const int numLevels,
            const int minSpatialLevel,
            const int minTemporalLevel,
            const bool loadInitMesh=false);

    //! Default constructor
    ~ Solver ();

    void setConfigParams();

    void setMeshHierarchy();

    void setDiscretisation();

    void setObserver();

    //! Runs the solver and computes the solution error
    std::tuple<int, double, double, Eigen::VectorXd>
    operator()();

    //! Initializes, assembles and then solves
    void run();

    void initialize();

    void assembleSystem();
    void assembleRhs();

    void solve();
    void solve(const Vector&, Vector&);
    void setPardisoSolver();

    //! Computes error
    Eigen::VectorXd evalError();
    Eigen::VectorXd evalError(BlockVector&);

private:
    double getMeshwidthOfFinestSpatialMesh();

    int getNumDofs() {
        return m_solutionHandler->getDataSize().Sum();
    }

public:
    std::shared_ptr<sparseHeat::LsqSparseXtFem> getDiscretisation(){
        return m_disc;
    }

    /*
    //! Releases allocated memory
    void finalize () const;

    //! Computes the solution error
    Eigen::VectorXd computeError(BlockVector *W) const;

    //! Evaluates the temperature solution at end time
    std::shared_ptr <GridFunction>
    getTemperatureSolAtEndTime (BlockVector *W) const;

    //! Evaluates the flux solution at end time
    std::shared_ptr <GridFunction>
    getFluxSolAtEndTime (BlockVector *W) const;

    //! Sets the temperature solution at end time
    void setTemperatureSolAtEndTime ();

    //! Sets the flux solution at end time
    void setFluxSolAtEndTime ();*/

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

    bool m_boolError;
    std::string m_errorType;

    double m_endTime;

    std::shared_ptr<mymfem::NestedMeshHierarchy> m_spatialMeshHierarchy;

    std::shared_ptr<sparseHeat::LsqSparseXtFem> m_disc;
//    std::shared_ptr<heat::Observer> m_observer;

    std::unique_ptr<sparseHeat::SolutionHandler> m_solutionHandler;
    std::shared_ptr <BlockVector> m_rhs;

    BlockOperator *m_systemOp = nullptr;
    SparseMatrix *m_systemMat = nullptr;

    std::shared_ptr <GridFunction> m_temperature;
    std::shared_ptr <GridFunction> m_flux;

#ifdef PARDISO_HPP
    std::unique_ptr<PardisoSolver> m_pardisoSolver;
#endif
};

}

#endif // SPARSE_HEAT_SOLVER_HPP
