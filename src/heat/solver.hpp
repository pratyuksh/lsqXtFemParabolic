#ifndef HEAT_SOLVER_HPP
#define HEAT_SOLVER_HPP

#include "mfem.hpp"

#include "../core/config.hpp"
#include "../pardiso/pardiso.hpp"

#include "../mymfem/utilities.hpp"

#include "test_cases_factory.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "solution_handler.hpp"


namespace heat{

/**
 * @brief Space-time solver for the heat equation
 *
 * Uses the least-squares FEM.
 */
class Solver
{
public:
    Solver (const nlohmann::json& config,
            std::shared_ptr<heat::TestCases>& testCase,
            std::string meshDir, int spatialLevel, int temporalLevel,
            bool loadInitMesh=false);

    Solver (const nlohmann::json& config,
            std::shared_ptr<heat::TestCases>& testCase,
            std::string meshDir, int spatialLevel,
            bool loadInitMesh=false);

    void setConfigParams();

    void setMeshes();

    void setDiscretisation();

    ~ Solver ();

    //! Sets perturbation in the test case
    void setPerturbation(double);

    void run();
    std::pair<double, int> runAndMeasurePerformanceMetrics();

    void initialize ();

    void assembleSystem();
    void assembleRhs();

    void solve ();
    std::pair<double, int> solveAndMeasurePerformanceMetrics ();
    std::pair<double, int> solve (const mfem::Vector&, mfem::Vector&);
    void setPardisoSolver();
    void finalizePardisoSolver();

    std::shared_ptr<mfem::Mesh> getTemporalMesh() const {
        return m_temporalMesh;
    }

    std::shared_ptr<mfem::Mesh> getSpatialMesh() const {
        return m_spatialMesh;
    }

    std::shared_ptr<heat::TestCases> getTestCase() const {
        return m_testCase;
    }

    std::shared_ptr<heat::LsqXtFem> getDiscretisation() const {
        return m_disc;
    }

    std::shared_ptr<heat::SolutionHandler> getSolutionHandler() const {
        return m_solutionHandler;
    }

    std::pair<double, double> getMeshwidths() const
    {
        double ddum;

        double htMin, htMax;
        m_temporalMesh->GetCharacteristics(htMin, htMax, ddum, ddum);

        double hxMin, hxMax;
        m_spatialMesh->GetCharacteristics(hxMin, hxMax, ddum, ddum);

        return {std::move(htMax), std::move(hxMax)};
    }

    int getNumDofs() const {
        return m_solutionHandler->getDataSize().Sum();
    }

private:
    const nlohmann::json& m_config;

    std::shared_ptr<heat::TestCases> m_testCase;

    std::string m_meshDir;
    std::string m_meshElemType;

    int m_spatialLevel, m_initSpatialLevel;
    int m_temporalLevel;
    bool m_loadInitMesh;

    std::string m_discType;
    std::string m_linearSolver;
    
    double m_endTime;

    std::shared_ptr<mfem::Mesh> m_temporalMesh;
    std::shared_ptr<mfem::Mesh> m_spatialMesh;

    std::shared_ptr<heat::LsqXtFem> m_disc;

    bool m_boolError;
    std::string m_errorType;
    
    mfem::Array<int> m_blockOffsets;

    mfem::BlockOperator *m_systemOp = nullptr;
    mfem::SparseMatrix *m_systemMat = nullptr;

    std::shared_ptr<heat::SolutionHandler> m_solutionHandler;
    std::unique_ptr <mfem::BlockVector> m_rhs;

#ifdef PARDISO_HPP
    std::unique_ptr<PardisoSolver> m_pardisoSolver;
#endif
};

}

#endif // HEAT_SOLVER_HPP
