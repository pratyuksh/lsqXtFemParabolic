#ifndef HEAT_SOLVER_HPP
#define HEAT_SOLVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include <Eigen/Core>

#include "../core/config.hpp"
#include "../pardiso/pardiso.hpp"

#include "../mymfem/utilities.hpp"

#include "test_cases_factory.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "observer.hpp"


namespace heat{

/**
 * @brief Space-time solver for the heat equation
 *
 * Uses the least-squares FEM.
 */
class Solver
{
public:
    //! Constructor with config JSON
    Solver (const nlohmann::json& config);

    //! Constructor with config JSON and test case
    Solver (const nlohmann::json& config,
            std::shared_ptr<heat::TestCases>& testCase);

    //! Constructor with config JSON and space-time mesh info
    Solver (const nlohmann::json& config,
            std::string mesh_dir, const int lx, const int lt,
            const bool load_init_mesh=false);

    //! Constructor with config JSON, test case and space-time mesh info
    Solver (const nlohmann::json& config,
            std::shared_ptr<heat::TestCases>& testCase,
            std::string mesh_dir, const int lx, const int lt,
            const bool load_init_mesh=false);

    //! Constructor with config JSON and spatial mesh info
    Solver (const nlohmann::json& config,
            std::string mesh_dir, const int lx,
            const bool load_init_mesh=false);

    //! Default constructor
    ~ Solver ();

    //! Sets auxiliary variables, meshes, discretization and observer
    void set(std::string mesh_dir, const int lx, const int lt=-2,
             const bool load_init_mesh=false);

    //! Sets mesh related auxiliary variables
    void set(std::string mesh_dir,
             const bool load_init_mesh=false);
    //! Sets meshes, discretization and observer
    void set(int lx, int lt=-2);
    //! Sets meshes
    void setMeshes(int lx, int lt=-2);

    //! Sets perturbation in the test case
    void setPerturbation(double);

    //! Initializes and then runs the solver,
    //! used for the full-tensor version
    std::tuple<int, double, double, Eigen::VectorXd>
    operator()();

    // Initializes and then runs the solver,
    // used for the sparse-grids handler
    // void operator()(std::unique_ptr<BlockVector>&);

    //! Runs the solver, assumes the solver has already been initialized
    void run();

    double computeObs();
    double computeQoi();
    double computeXtQoi();
    double getObs();
    double getQoi();
    double getXtQoi();
    
    //! Initializes the solver
    void init ();

private:
    //! Assembles the right-hand side
    void assembleRhs();
    //! Assembles the system
    void assembleSystem();

public:
    //! Solves the linear system resulting from the discretisation
    void solve (BlockVector *W, BlockVector* B) const;
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
    void setFluxSolAtEndTime ();

    //! Initializes and computes the projection,
    //! Used for the full-tensor version
    std::tuple<int, double, double, Eigen::VectorXd>
    projection();

    // Initializes and computes the projection,
    // Used for the sparse-tensor version
    // void projection(std::unique_ptr<BlockVector>&);

    //! Initializes projection
    void initProjection ();
    //! Assembles projection right-hand side
    void assembleProjectionRhs();
    //! Assembles the system matrix for projection
    void assembleProjectionSystem();

    //! Returns the temporal mesh
    inline std::shared_ptr<Mesh> getTMesh() const {
        return m_tMesh;
    }

    //! Returns the spatial mesh
    inline std::shared_ptr<Mesh> getXMesh() const {
        return m_xMesh;
    }

    //! Returns the test case
    inline std::shared_ptr<heat::TestCases>
    getTestCase() const {
        return m_testCase;
    }

    //! Returns the discretisation
    inline std::shared_ptr<heat::LsqXtFEM>
    getDiscretisation() const {
        return m_discr;
    }

    //! Returns the observer
    inline std::shared_ptr<heat::Observer>
    getObserver() const {
        return m_observer;
    }

    //! Returns the mesh characteristics
    inline std::pair<double, double>
    getMeshChars() const
    {
        double htMin, htMax, kappatMin, kappatMax;
        m_tMesh->GetCharacteristics(htMin, htMax, kappatMin, kappatMax);

        double hxMin, hxMax, kappaxMin, kappaxMax;
        m_xMesh->GetCharacteristics(hxMin, hxMax, kappaxMin, kappaxMax);

        return {std::move(htMax), std::move(hxMax)};
    }

    //! Returns the total number of degrees of freedom in space-time
    inline int getNdofs() const {
        return std::move(m_W->Size());
    }

private:
    const nlohmann::json& m_config;
    std::string m_meshDir;
    std::string m_meshElemType;
    std::string m_xDiscrType;
    std::string m_linearSolver;
    
    double m_endTime;
    std::shared_ptr<Mesh> m_tMesh;
    std::shared_ptr<Mesh> m_xMesh;
    bool m_loadInitMesh;
    int m_lx0;
    
    std::shared_ptr<heat::TestCases> m_testCase;
    std::shared_ptr<heat::LsqXtFEM> m_discr;
    std::shared_ptr<heat::Observer> m_observer;

    bool m_boolError;
    std::string m_errorType;
    
    Array<int> m_blockOffsets;

    BlockOperator *m_heatOp = nullptr;
    SparseMatrix *m_heatMat = nullptr;
    SparseMatrix *m_heatProjMat = nullptr;

    std::unique_ptr <BlockVector> m_W;
    std::unique_ptr <BlockVector> m_B;
    bool m_firstPass = true;

    std::shared_ptr <GridFunction> m_temperature;
    std::shared_ptr <GridFunction> m_flux;

#ifdef PARDISO_HPP
    std::unique_ptr<PardisoSolver> m_pardisoSolver;
#endif
    mutable bool m_pardisoFinalized = true;
};

}

#endif // HEAT_SOLVER_HPP
