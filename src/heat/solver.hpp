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

class Solver
{
public:
    Solver (const nlohmann::json& config);

    Solver (const nlohmann::json& config,
                std::shared_ptr<heat::TestCases>& testCase);

    Solver (const nlohmann::json& config,
                std::string mesh_dir,
                const int lx,
                const int lt,
                const bool load_init_mesh=false);

    Solver (const nlohmann::json& config,
                std::shared_ptr<heat::TestCases>& testCase,
                std::string mesh_dir,
                const int lx,
                const int lt,
                const bool load_init_mesh=false);

    Solver (const nlohmann::json& config,
                std::string mesh_dir,
                const int lx,
                const bool load_init_mesh=false);

    ~ Solver ();

    // sets meshes, discretization, observer
    void set(std::string mesh_dir,
             const int lx,
             const int lt=-2,
             const bool load_init_mesh=false);

    void set(std::string mesh_dir,
             const bool load_init_mesh=false);
    void set(int lx, int lt=-2);
    void setMeshes(int lx, int lt=-2);

    // sets perturbation in the test case
    void setPerturbation(double);

    // used for the full-tensor version
    std::tuple<int, double, double, Eigen::VectorXd>
    operator()();

    // used for the sparse-grids handler
    void operator()(std::unique_ptr<BlockVector>&);

    void run();
    double computeObs();
    double computeQoi();
    double computeXtQoi();
    double getObs();
    double getQoi();
    double getXtQoi();
    
    void init ();

private:
    void assembleRhs();
    void assembleSystem();

public:
    void solve (BlockVector *W, BlockVector* B) const;
    void finalize () const;

    // error computation
    Eigen::VectorXd computeError(BlockVector *W) const;

    // solution at end-time
    std::shared_ptr <GridFunction>
    getTemperatureSolAtEndTime (BlockVector *W) const;

    std::shared_ptr <GridFunction>
    getFluxSolAtEndTime (BlockVector *W) const;

    void setTemperatureSolAtEndTime ();

    void setFluxSolAtEndTime ();

    // for L2-projection
    void projection(std::unique_ptr<BlockVector>&);
    std::tuple<int, double, double, Eigen::VectorXd>
    projection();

    void initProjection ();
    void assembleProjectionRhs();
    void assembleProjectionSystem();

    inline std::shared_ptr<Mesh> getTMesh() const {
        return m_tMesh;
    }

    inline std::shared_ptr<Mesh> getXMesh() const {
        return m_xMesh;
    }

    inline std::shared_ptr<heat::TestCases>
    getTestCase() const {
        return m_testCase;
    }

    inline std::shared_ptr<heat::LsqXtFEM>
    getDiscretisation() const {
        return m_discr;
    }

    inline std::shared_ptr<heat::Observer>
    getObserver() const {
        return m_observer;
    }

    inline std::pair<double, double>
    getMeshChars() const
    {
        double ht_min, ht_max, kappat_min, kappat_max;
        m_tMesh->GetCharacteristics(ht_min, ht_max,
                                    kappat_min, kappat_max);

        double hx_min, hx_max, kappax_min, kappax_max;
        m_xMesh->GetCharacteristics(hx_min, hx_max,
                                    kappax_min, kappax_max);

        return {std::move(ht_max), std::move(hx_max)};
    }

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
