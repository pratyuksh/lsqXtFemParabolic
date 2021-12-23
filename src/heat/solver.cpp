#include "solver.hpp"
#include "utilities.hpp"

#include <iostream>
#include <chrono>

using namespace mfem;


// Constructor with config JSON, test case and space-time mesh info
heat::Solver
:: Solver (const nlohmann::json& config,
           std::shared_ptr<heat::TestCases>& testCase,
           std::string meshDir,
           int spatialLevel,
           int temporalLevel,
           bool loadInitMesh)
    : m_config (config),
      m_testCase(testCase),
      m_meshDir (meshDir),
      m_spatialLevel (spatialLevel),
      m_temporalLevel (temporalLevel),
      m_loadInitMesh (loadInitMesh)
{
    setConfigParams();
    setMeshes();
    setDiscretisation();
}

// Constructor with config JSON, test case and spatial mesh info
heat::Solver
:: Solver (const nlohmann::json& config,
           std::shared_ptr<heat::TestCases>& testCase,
           std::string meshDir,
           int spatialLevel,
           bool loadInitMesh)
    : Solver(config, testCase, meshDir, spatialLevel, -2, loadInitMesh) {}

heat::Solver :: ~ Solver () {}

void heat::Solver
:: setConfigParams()
{
    m_endTime = m_config["end_time"];

    m_meshElemType = "tri";
    if (m_config.contains("mesh_elem_type")) {
        m_meshElemType = m_config["mesh_elem_type"];
    }

    m_initSpatialLevel = 0;
    if (m_loadInitMesh) {
        if (m_config.contains("init_mesh_level")) {
            m_initSpatialLevel = m_config["init_mesh_level"];
        }
    }

    m_discType = "H1Hdiv";
    if (m_config.contains("discretisation_type")) {
        m_discType = m_config["discretisation_type"];
    }

    m_linearSolver = "pardiso";
    if (m_config.contains("linear_solver")) {
        m_linearSolver = m_config["linear_solver"];
    }

    m_boolError = false;
    if (m_config.contains("eval_error")) {
        m_boolError = m_config["eval_error"];
    }

    m_errorType = "natural";
    if (m_config.contains("error_type")) {
        m_errorType = m_config["error_type"];
    }
}

void heat::Solver
:: setMeshes()
{
    // mesh in space
    if (m_loadInitMesh) { // load initial mesh and refine
        int numRefinements = m_spatialLevel - m_initSpatialLevel;
        const std::string meshFile
                        = m_meshDir+"/"+m_meshElemType
                +"_mesh_l"
                +std::to_string(m_initSpatialLevel)+".mesh";
#ifndef NDEBUG
        std::cout << "  Initial mesh file: " << meshFile << std::endl;
        std::cout << "  Number of uniform refinements: "
                  << numRefinements << std::endl;
#endif
        m_spatialMesh = std::make_shared<Mesh>(meshFile.c_str());
        for (int k=0; k<numRefinements; k++) {
            m_spatialMesh->UniformRefinement();
        }
    }
    else {
        const std::string meshFile =
                m_meshDir+"/mesh_l"+std::to_string(m_spatialLevel)+".mesh";
#ifndef NDEBUG
        std::cout << "  Mesh file: "
                  << meshFile << std::endl;
#endif
        m_spatialMesh = std::make_shared<Mesh>(meshFile.c_str());
    }

    // mesh in time
    if (m_temporalLevel == -2) { // set time-level acc. to spatial-level
        double hxMin, hxMax, ddum;
        m_spatialMesh->GetCharacteristics(hxMin, hxMax, ddum, ddum);
        m_temporalLevel = int(std::ceil(-std::log2(hxMax/m_endTime)));
    }
    int Nt = static_cast<int>(std::pow(2, m_temporalLevel));
    m_temporalMesh = std::make_shared<Mesh>(Nt, m_endTime);
}

void heat::Solver
:: setDiscretisation()
{
    if (m_discType == "H1Hdiv") {
        m_disc = std::make_shared<heat::LsqXtFemH1Hdiv>
                (m_config, m_testCase);
    }
    else if (m_discType == "H1H1") {
        m_disc = std::make_shared<heat::LsqXtFemH1H1>
                (m_config, m_testCase);
    }
    else {
        std::cout << "Unknown discretisation!" << std::endl;
        abort();
    }

    //m_discr->setFESpacesAndSpatialBoundaryDofs(m_tMesh, m_xMesh);
    m_disc->setFeSpacesAndBlockOffsetsAndSpatialBoundaryDofs
            (m_temporalMesh, m_spatialMesh);
    m_blockOffsets = m_disc->getBlockOffsets();
}

void heat::Solver
:: setPerturbation(double omega) {
    m_testCase->setPerturbation(omega);
}

void heat::Solver
:: run()
{
    initialize ();
    assembleSystem();
    assembleRhs();
    solve ();
}

std::pair<double, int> heat::Solver
:: runAndMeasurePerformanceMetrics()
{
    initialize ();
    assembleSystem();
    assembleRhs();
    return solveAndMeasurePerformanceMetrics();
}

// Initializes the solver
void heat::Solver
:: initialize ()
{
    // solution data handler
    auto temporalFeSpace = m_disc->getTemporalFeSpace();
    auto spatialFeSpaces = m_disc->getSpatialFeSpaces();
    m_solutionHandler
            = std::make_shared<heat::SolutionHandler>
            (temporalFeSpace, spatialFeSpaces);

    // variable to store rhs data
    m_rhs = std::make_unique<BlockVector>(m_blockOffsets);
}

void heat::Solver
:: assembleSystem()
{
//    m_discr->assembleSystem();
    m_disc->assembleSystemSubMatrices();

    if (m_linearSolver == "pardiso")
    {
//        m_discr->buildSystemMatrix();
        m_disc->buildUpperTriangleOfSystemMatrix();
        m_systemMat = m_disc->getSystemMatrix();
    }
    else if (m_linearSolver == "cg")
    {
        m_disc->buildSystemMatrix();
        m_systemMat = m_disc->getSystemMatrix();
    }
}

void heat::Solver
:: assembleRhs()
{
    (*m_rhs) = 0.0;
    m_disc->assembleRhs(m_rhs.get());
}

void heat::Solver
:: solve()
{
    // wrap rhs BlockVector as vector
    Vector wrapRhs(m_rhs->GetData(), m_rhs->Size());

    // wrap solution BlockVector as vector
    auto U = m_solutionHandler->getData();
    Vector wrapU(U->GetData(), U->Size());

    solve(wrapRhs, wrapU);
}

std::pair<double, int> heat::Solver
:: solveAndMeasurePerformanceMetrics()
{
    // wrap rhs BlockVector as vector
    Vector wrapRhs(m_rhs->GetData(), m_rhs->Size());

    // wrap solution BlockVector as vector
    auto U = m_solutionHandler->getData();
    Vector wrapU(U->GetData(), U->Size());

    return solve(wrapRhs, wrapU);
}

std::pair<double, int> heat::Solver
:: solve (const mfem::Vector& rhs, mfem::Vector& u)
{
    int memUsage = 0;

    auto start = std::chrono::high_resolution_clock::now();
    if (m_linearSolver == "pardiso")
    {
        setPardisoSolver();
        m_pardisoSolver->initialize(m_systemMat->Size(),
                                    m_systemMat->GetI(),
                                    m_systemMat->GetJ(),
                                    m_systemMat->GetData());
        m_pardisoSolver->factorize();
        m_pardisoSolver->solve(rhs.GetData(), u.GetData());
        memUsage = m_pardisoSolver->getMemoryUsage();
        finalizePardisoSolver();
    }
    else if (m_linearSolver == "cg")
    {
        int verbose = m_config["cg_verbose"];
        int maxIter = m_config["cg_maxIter"];
        double relTol = m_config["cg_relTol"];
        double absTol = m_config["cg_absTol"];

        auto M = new GSSmoother(*m_systemMat);
        u = 0.;
        PCG(*m_systemMat, *M, rhs, u, verbose, maxIter, relTol, absTol);
        delete M;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast
            <std::chrono::milliseconds>(end - start);
    double elapsedTime
            = (static_cast<double>(duration.count()))/1000;

    return {elapsedTime, memUsage};
}

void heat::Solver
:: setPardisoSolver()
{
//    int mtype = 1; // real structurally symmetric
        int mtype = 2; // real symmetric positive definite
//        int mtype = -2; // real symmetric indefinite
    m_pardisoSolver
            = std::make_unique<PardisoSolver>(mtype);
}

void heat::Solver
:: finalizePardisoSolver() {
    m_pardisoSolver->finalize();
}

// End of file
