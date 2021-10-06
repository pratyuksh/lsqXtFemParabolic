#include "solver.hpp"
#include "utilities.hpp"

#include <iostream>
#include <chrono>


// Constructor
heat::Solver :: Solver (const nlohmann::json& config)
    : m_config(config)
{
    m_endTime = config["end_time"];

    m_boolError = false;
    if (config.contains("eval_error")) {
        m_boolError = config["eval_error"];
    }

    m_errorType = "H1";
    if (config.contains("error_type")) {
        m_errorType = config["error_type"];
    }

    m_xDiscrType = "H1Hdiv";
    if (config.contains("discretisation_type")) {
        m_xDiscrType = config["discretisation_type"];
    }

    m_linearSolver = "pardiso";
    if (config.contains("linear_solver")) {
        m_linearSolver = config["linear_solver"];
    }
}

heat::Solver :: Solver (const nlohmann::json& config,
                          std::shared_ptr<heat::TestCases>&
                          testCase)
    : Solver (config)
{
    m_testCase = testCase;
}

heat::Solver :: Solver (const nlohmann::json& config,
                          std::string mesh_dir,
                          const int lx,
                          const int lt,
                          const bool load_init_mesh)
    : Solver (config)
{
    set(mesh_dir, lx, lt, load_init_mesh);
}

heat::Solver :: Solver (const nlohmann::json& config,
                          std::shared_ptr<heat::TestCases>&
                          testCase,
                          std::string mesh_dir,
                          const int lx,
                          const int lt,
                          const bool load_init_mesh)
    : Solver (config)
{
    m_testCase = testCase;
    set(mesh_dir, lx, lt, load_init_mesh);
}

heat::Solver :: Solver (const nlohmann::json& config,
                          std::string mesh_dir,
                          const int lx,
                          const bool load_init_mesh)
    : Solver (config)
{
    set(mesh_dir, lx, -2, load_init_mesh);
}

// Destructor
heat::Solver :: ~ Solver ()
{
    finalize();
}

// Sets test case, meshes, discretization, observer
void heat::Solver :: set (std::string mesh_dir,
                        const int lx,
                        const int lt,
                        const bool load_init_mesh)
{
    set(mesh_dir, load_init_mesh);
#ifdef MYVERBOSE
    std::cout << "\nRun Heat solver ..." << std::endl;
#endif
    set(lx, lt);
}

void heat::Solver :: set (std::string mesh_dir,
                        const bool load_init_mesh)
{
    m_loadInitMesh = load_init_mesh;
    m_meshDir = mesh_dir;

    m_lx0 = 0;
    if (m_loadInitMesh) {
        if (m_config.contains("init_mesh_level")) {
            m_lx0 = m_config["init_mesh_level"];
        }
    }

    m_meshElemType = "tri";
    if (m_config.contains("mesh_elem_type")) {
        m_meshElemType = m_config["mesh_elem_type"];
    }
}

void heat::Solver :: set(int lx, int lt)
{
    m_firstPass = true;

    // test case
    if (!m_testCase) {
        m_testCase = heat::makeTestCase(m_config);
    }

    // meshes
    setMeshes(lx, lt);

    // space-time discretisation
    m_discr.reset();
    if (m_xDiscrType == "H1Hdiv") {
        m_discr = std::make_shared<heat::LsqXtFemH1Hdiv>
                (m_config, m_testCase);
    }
    else if (m_xDiscrType == "H1H1") {
        m_discr = std::make_shared<heat::LsqXtFemH1H1>
                (m_config, m_testCase);
    }
    else {
        std::cout << "Unknown discretisation!" << std::endl;
        abort();
    }

    m_discr->set(m_tMesh, m_xMesh);
    m_blockOffsets = m_discr->getBlockOffsets();

    // observer
    m_observer.reset();
    m_observer = std::make_shared<heat::Observer>
            (m_config, lx);
    m_observer->set(m_testCase, m_discr);
}

// Sets meshes
void heat::Solver :: setMeshes(int lx, int lt)
{
    // mesh in space
    m_xMesh.reset();
    if (m_loadInitMesh) { // load initial mesh and refine
        int num_refinements = lx - m_lx0;
        const std::string mesh_file
                        = m_meshDir+"/"+m_meshElemType
                +"_mesh_l"
                +std::to_string(m_lx0)+".mesh";
#ifdef MYVERBOSE
        std::cout << "  Initial mesh file: " << mesh_file << std::endl;
        std::cout << "  Number of uniform refinements: "
                  << num_refinements << std::endl;
#endif
        m_xMesh = std::make_shared<Mesh>(mesh_file.c_str());
        for (int k=0; k<num_refinements; k++) {
            m_xMesh->UniformRefinement();
        }
    }
    else {
        const std::string mesh_file =
                m_meshDir+"/mesh_l"+std::to_string(lx)+".mesh";
#ifdef MYVERBOSE
        std::cout << "  Mesh file: "
                  << mesh_file << std::endl;
#endif
        m_xMesh = std::make_shared<Mesh>(mesh_file.c_str());
    }

    // mesh in time
    m_tMesh.reset();
    if (lt == -2) {
        double hx_min, hx_max, kappax_min, kappax_max;
        m_xMesh->GetCharacteristics(hx_min, hx_max,
                                    kappax_min, kappax_max);

        // set time-level acc. to spatial-level
        lt = int(std::ceil(-std::log2(hx_max/m_endTime)));
    }
    int Nt = static_cast<int>(std::pow(2,lt));
    m_tMesh = std::make_shared<Mesh>(Nt, m_endTime);
}

// Sets perturbation in the test case
void heat::Solver :: setPerturbation(double omega)
{
    m_testCase->setPerturbation(omega);
}

// Runs the solver; used for the full-tensor version
std::tuple<int, double, double, Eigen::VectorXd> heat::Solver
:: operator() ()
{
    // run solver
    run();

    // visualization at end time
    auto u = getTemperatureSolAtEndTime(m_W.get());
    (*m_observer)(u);

    // compute error
    auto errSol = computeError(m_W.get());

    // max meshwidth in space and time
    double ht_max, hx_max=-1;
    std::tie(ht_max, hx_max) = getMeshChars();

    return {std::move(getNdofs()),
                std::move(ht_max), std::move(hx_max),
                std::move(errSol) };
}

// Runs the solver; used for the sparse-grids handler
void heat::Solver
:: operator()(std::unique_ptr<BlockVector>& W)
{
    // MFEM solution variables
    W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM rhs variables
    m_B = std::make_unique<BlockVector>(m_blockOffsets);

    // initialize
    init ();

    // solve
    solve (W.get(), m_B.get());

    // visualization at end time
    auto u = getTemperatureSolAtEndTime(W.get());
    (*m_observer)(u);
}

// Runs the solver
void heat::Solver
:: run()
{
    if (m_firstPass) {
        // MFEM solution variables
        m_W.reset();
        m_W = std::make_unique<BlockVector>(m_blockOffsets);

        // MFEM rhs variables
        m_B.reset();
        m_B = std::make_unique<BlockVector>(m_blockOffsets);
    }

    // initialize
    init ();

    // solve
    solve (m_W.get(), m_B.get());
}

// Runs the solver
// and evaluates the flux observation functional
double heat::Solver
:: computeObs()
{
    run();
    //auto u = get_temperatureSol_at_endTime(m_W.get());
    //(*m_observer)(u);
    setTemperatureSolAtEndTime();
    //set_fluxSol_at_endTime();
    return std::move(getObs());
}

// Runs the solver
// and evaluates the flux quantity of interest
double heat::Solver
:: computeQoi()
{
    run();
    //auto u = get_temperatureSol_at_endTime(m_W.get());
    //(*m_observer)(u);
    setTemperatureSolAtEndTime();
    //set_fluxSol_at_endTime();
    return std::move(getQoi());
}

double heat::Solver
:: computeXtQoi()
{
    run();
    return std::move(getXtQoi());
}

// Evaluates the flux observation functional
double heat::Solver
:: getObs()
{
    auto fluxObs = m_observer->evalObs(m_temperature);
    return std::move(fluxObs);
}

// Evaluates the quantity of interest functional
double heat::Solver
:: getQoi()
{
    auto fluxQoi = m_observer->evalQoi(m_temperature);
    return std::move(fluxQoi);
}

double heat::Solver
:: getXtQoi()
{
    auto fluxQoi = m_observer->evalXtQoi(m_W->GetBlock(0));
    return std::move(fluxQoi);
}

// Initializes the solver
void heat::Solver :: init ()
{
    assembleRhs();
    assembleSystem();
}

// Assembles the right-hand side
void heat::Solver :: assembleRhs()
{
    if (m_firstPass) {
        (*m_B) = 0.0;
        m_discr->assembleRhs(m_B.get());
        m_firstPass = false;
    }
}

// Assembles the system
void heat::Solver :: assembleSystem()
{
    m_discr->assembleSystem();

    if (!m_pardisoFinalized) {
#ifdef MYVERBOSE
        std::cout << "\nRelease old Pardiso memory"
                  << std::endl;
#endif
        m_pardisoSolver->finalize();
        m_pardisoFinalized = true;
    }

    if (m_linearSolver == "pardiso")
    {
        //m_discr->build_system_matrix();
        m_discr->buildSystemMatrixUpTr();
        m_heatMat = m_discr->getHeatMat();
        //auto start = std::chrono::high_resolution_clock::now();
        // set Pardiso solver
        //int mtype = 1; // real structurally symmetric
        //int mtype = 2; // real symmetric positive definite
        int mtype = -2; // real symmetric indefinite
        m_pardisoSolver.reset();
        m_pardisoSolver = std::make_unique<PardisoSolver>(mtype);
        m_pardisoSolver->initialize(m_heatMat->Size(),
                                    m_heatMat->GetI(),
                                    m_heatMat->GetJ(),
                                    m_heatMat->GetData());
        m_pardisoSolver->factorize();
        m_pardisoFinalized = false;
        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast
                <std::chrono::milliseconds>(end - start);
        double elapsedTime
                = (static_cast<double>(duration.count()))/1000;
        std::cout << "Factorization time: " << elapsedTime << std::endl;
        */
    }
    else if (m_linearSolver == "cg")
    {
        m_discr->buildSystemMatrix();
        m_heatMat = m_discr->getHeatMat();
    }
}

// Solves the linear system
void heat::Solver :: solve (BlockVector *W,
                          BlockVector *B) const
{
    if (m_linearSolver == "pardiso")
    {
        //auto start = std::chrono::high_resolution_clock::now();
        m_pardisoSolver->solve(B->GetData(), W->GetData());
        /*auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast
                <std::chrono::milliseconds>(end - start);
        double elapsedTime
                = (static_cast<double>(duration.count()))/1000;
        std::cout << "Solve time: " << elapsedTime << std::endl;*/
    }
    else if (m_linearSolver == "cg")
    {
        int verbose = m_config["cg_verbose"];
        int maxIter = m_config["cg_maxIter"];
        double relTol = m_config["cg_relTol"];
        double absTol = m_config["cg_absTol"];

        GSSmoother *M = new GSSmoother(*m_heatMat);
        (*W) = 0;
        PCG(*m_heatMat, *M, *B, *W, verbose, maxIter, relTol, absTol);
        delete M;
    }
}

// Releases memory
void heat::Solver :: finalize() const
{
    if (!m_pardisoFinalized) {
#ifdef MYVERBOSE
        std::cout << "\nRelease Pardiso Memory" << std::endl;
#endif
        m_pardisoSolver->finalize();
        m_pardisoFinalized = true;
    }
}

// Computes the numerical solution error
Eigen::VectorXd heat::Solver
:: computeError(BlockVector *W) const
{
    Eigen::VectorXd errSol(2);
    errSol.setZero();

    if (m_boolError)
    {
        if (m_errorType=="H1") {
#ifdef MYVERBOSE
            std::cout << "  Computing H1 error "
                         "at end time..."
                      << std::endl;
#endif
            auto u = getTemperatureSolAtEndTime(W);

            double erruL2, uEL2, erruH1, uEH1;
            std::tie (erruL2, uEL2, erruH1, uEH1)
                    = m_observer->evalXH1Error(u, m_endTime);
#ifdef MYVERBOSE
            std::cout << "\ttemperature L2 error: "
                      << erruL2 << "\t" << uEL2 << std::endl;
            std::cout << "\ttemperature H1 error: "
                      << erruH1 << "\t" << uEH1 << std::endl;
#endif
            errSol(0) = erruL2/uEL2;
            errSol(1) = erruH1/uEH1;
        }
        else if (m_errorType == "L2H1") {
#ifdef MYVERBOSE
            std::cout << "  Computing L2H1 "
                         "space-time error..."
                      << std::endl;
#endif
            double erruL2L2, uEL2L2, erruL2H1, uEL2H1;
            std::tie (erruL2L2, uEL2L2, erruL2H1, uEL2H1)
                    = m_observer->evalXtL2H1Error
                    (W->GetBlock(0));
#ifdef MYVERBOSE
            std::cout << "\ttemperature L2L2 error: "
                      << erruL2L2 << "\t"
                      << uEL2L2 << std::endl;
            std::cout << "\ttemperature L2H1 error: "
                 << erruL2H1 << "\t"
                 << uEL2H1 << std::endl;
#endif
            errSol(0) = erruL2L2/uEL2L2;
            errSol(1) = erruL2H1/uEL2H1;
        }
        else if (m_errorType == "lsq") {
//#ifdef MYVERBOSE
            std::cout << "  Computing space-time "
                         "least-squares error..."
                      << std::endl;
//#endif
            double errPde, errFlux, errIc;
            std::tie (errPde, errFlux, errIc)
                    = m_observer->evalXtLsqError
                    (W->GetBlock(0), W->GetBlock(1));
//#ifdef MYVERBOSE
            std::cout << "\tLsq Pde error: "
                      << errPde << std::endl;
            std::cout << "\tLsq Flux error: "
                      << errFlux << std::endl;
            std::cout << "\tLsq Ic error: "
                      << errIc << std::endl;
//#endif
            errSol.resize(3);
            errSol(0) = errPde;
            errSol(1) = errFlux;
            errSol(2) = errIc;
        }
        else if (m_errorType == "natural") {
//#ifdef MYVERBOSE
            std::cout << "  Computing error "
                         "in natural norm..."
                      << std::endl;
//#endif
            double erruL2H1, errqL2L2, errUDiv;
            std::tie (erruL2H1, errqL2L2, errUDiv)
                    = m_observer->evalXtError
                    (W->GetBlock(0), W->GetBlock(1));
//#ifdef MYVERBOSE
            std::cout << "\tu L2H1 error: "
                      << erruL2H1 << std::endl;
            std::cout << "\tq L2L2 error: "
                      << errqL2L2 << std::endl;
            std::cout << "\tDivergence error: "
                      << errUDiv << std::endl;
//#endif
            errSol.resize(3);
            errSol(0) = erruL2H1;
            errSol(1) = errqL2L2;
            errSol(2) = errUDiv;
        }
        else {
            std::cout << "Unknown error type!\n";
            abort();
        }
    }

    return errSol;
}

// Evaluates the temperature solution vector at end time T
std::shared_ptr <GridFunction> heat::Solver
:: getTemperatureSolAtEndTime(BlockVector *W) const
{
    // auxiliary variables
    FiniteElementSpace* tFespace
                = m_discr->getTFespace();
    Array<FiniteElementSpace*> xFespaces
            = m_discr->getXFespaces();

    // dofs for the time element (T-h, T)
    Array<int> tVdofs;
    const FiniteElement *tFe
            = tFespace->GetFE(m_tMesh->GetNE()-1);
    tFespace->GetElementVDofs(m_tMesh->GetNE()-1, tVdofs);

    Vector tShape(tVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    tFe->CalcShape(ip, tShape);

    // temperature solution at T
    Vector uSol;
    std::shared_ptr <GridFunction> u
            = std::make_shared<GridFunction>(xFespaces[0]);
    u->GetTrueDofs(uSol);
    buildXSolFG(W->GetBlock(0), tShape, tVdofs, uSol);

    return u;
}

// Evaluates the flux solution vector at end time T
std::shared_ptr <GridFunction> heat::Solver
:: getFluxSolAtEndTime(BlockVector *W) const
{
    // auxiliary variables for visualization
    // and error computation at time T
    FiniteElementSpace* tFespace
                = m_discr->getTFespace();
    Array<FiniteElementSpace*> xFespaces
            = m_discr->getXFespaces();

    // dofs for the time element (T-h, T)
    Array<int> tVdofs;
    const FiniteElement *tFe
            = tFespace->GetFE(m_tMesh->GetNE()-1);
    tFespace->GetElementVDofs(m_tMesh->GetNE()-1, tVdofs);

    Vector tShape(tVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    tFe->CalcShape(ip, tShape);

    // flux solution at T
    Vector qSol;
    std::shared_ptr <GridFunction> q
            = std::make_shared<GridFunction>(xFespaces[1]);
    q->GetTrueDofs(qSol);
    buildXSolFG(W->GetBlock(1), tShape, tVdofs, qSol);

    return q;
}

// Sets the temperature solution vector at end time T
void heat::Solver
:: setTemperatureSolAtEndTime()
{
    m_temperature = getTemperatureSolAtEndTime(m_W.get());
}

// Sets the flux solution vector at end time T
void heat::Solver
:: setFluxSolAtEndTime()
{
    m_flux = getFluxSolAtEndTime(m_W.get());
}

// Routines for computing L2-projection
std::tuple<int, double, double, Eigen::VectorXd> heat::Solver
:: projection ()
{
    // MFEM solution variables
    m_W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM rhs variables
    m_B = std::make_unique<BlockVector>(m_blockOffsets);

    // initialize_projection
    initProjection();

    // solve
    solve (m_W.get(), m_B.get());

    // visualization at end time
    auto u = getTemperatureSolAtEndTime(m_W.get());
    (*m_observer)(u);

    // compute error
    auto errSol = computeError(m_W.get());

    // max meshwidth in space and time
    double ht_max, hx_max;
    std::tie(ht_max, hx_max) = getMeshChars();

    return {std::move(getNdofs()),
                std::move(ht_max), std::move(hx_max),
                std::move(errSol) };
}

void heat::Solver
:: projection (std::unique_ptr<BlockVector>& W)
{
    // MFEM solution variables
    W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM rhs variables
    m_B = std::make_unique<BlockVector>(m_blockOffsets);

    // initialize projection
    initProjection();

    // compute projection
    solve (W.get(), m_B.get());
}

void heat::Solver :: initProjection()
{
    assembleProjectionRhs();
    assembleProjectionSystem();
}

void heat::Solver :: assembleProjectionRhs()
{
    if (m_firstPass) {
        (*m_B) = 0.0;
        m_discr->assembleProjectionRhs(m_B.get());
        m_firstPass = false;
    }
}

void heat::Solver :: assembleProjectionSystem()
{
    m_discr->assembleProjector();
    m_discr->buildProjectorMatrix();
    m_heatProjMat = m_discr->getHeatProjectionMat();

    // set Pardiso solver
    int mtype = 2; // real symmetric positive definite
    m_pardisoSolver = std::make_unique<PardisoSolver>(mtype);
    m_pardisoSolver->initialize(m_heatProjMat->Size(),
                                m_heatProjMat->GetI(),
                                m_heatProjMat->GetJ(),
                                m_heatProjMat->GetData());
    m_pardisoSolver->factorize();
    m_pardisoFinalized = false;
}

// End of file
