#include "solver.hpp"
#include "utilities.hpp"

#include <iostream>
#include <chrono>

using namespace mfem;


sparseHeat::Solver
:: Solver (const nlohmann::json& config,
           std::shared_ptr<heat::TestCases>& testCase,
           std::string meshDir,
           int numLevels,
           int minSpatialLevel,
           int minTemporalLevel,
           bool loadInitMesh)
    : m_config (config),
      m_testCase(testCase),
      m_meshDir (meshDir),
      m_numLevels (numLevels),
      m_minSpatialLevel (minSpatialLevel),
      m_minTemporalLevel (minTemporalLevel),
      m_loadInitMesh (loadInitMesh)
{
    m_maxTemporalLevel = m_minTemporalLevel + m_numLevels - 1;

    setConfigParams();
    setMeshHierarchy();
    setDiscretisation();
}

sparseHeat::Solver
:: ~ Solver ()
{}

void sparseHeat::Solver
:: setConfigParams()
{
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(m_config, "end_time",
                                        m_endTime, 1);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(m_config, "mesh_elem_type",
                                        m_meshElemType, "tri");

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(m_config, "discretisation_type",
                                        m_discType, "H1Hdiv");

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(m_config, "linear_solver",
                                        m_linearSolver, "pardiso");
}

void sparseHeat::Solver
:: setMeshHierarchy()
{
    m_spatialMeshHierarchy
            = std::make_shared<mymfem::NestedMeshHierarchy>();

    if (m_loadInitMesh) // load initial mesh and refine
    {
        const std::string meshFile
                = m_meshDir+"/"+m_meshElemType
                +"_mesh_l0.mesh";
#ifndef NDEBUG
        std::cout << "  Initial mesh file: "
                  << meshFile << std::endl;
        std::cout << "  Minimum refinement level: "
                  << m_minSpatialLevel << std::endl;
        std::cout << "  Maximum refinement level: "
                  << m_minSpatialLevel + m_numLevels-1 << std::endl;
#endif
        // add coarsest spatial mesh to hierarchy
        auto xMesh = std::make_shared<Mesh>(meshFile.c_str());
        for (int k=0; k<m_minSpatialLevel; k++) {
            xMesh->UniformRefinement();
        }
        m_spatialMeshHierarchy->addMesh(xMesh);

        // refine and add to spatial mesh hierarchy
        auto xMeshOld = xMesh;
        for (int k=1; k<m_numLevels; k++) {
            auto xMeshNew = std::make_shared<Mesh>(*xMeshOld);
            xMeshNew->UniformRefinement();
            m_spatialMeshHierarchy->addMesh(xMeshNew);
            xMeshOld = xMeshNew;
        }
    }
    else {
        for (int k=0; k<m_numLevels; k++)
        {
            const std::string meshFile
                    = m_meshDir+"/mesh_l"
                    +std::to_string(m_minSpatialLevel + k)
                    +".mesh";
#ifndef NDEBUG
            std::cout << "  Mesh file: "
                      << meshFile << std::endl;
#endif
            // add current mesh to mesh hierarchy
            auto xMesh = std::make_shared<Mesh>(meshFile.c_str());
            m_spatialMeshHierarchy->addMesh(xMesh);
        }
    }
    m_spatialMeshHierarchy->finalize();
}

void sparseHeat::Solver
:: setDiscretisation()
{
    if (m_discType == "H1Hdiv") {
        m_disc = std::make_shared<sparseHeat::LsqSparseXtFemH1Hdiv>
                (m_config, m_testCase,
                 m_numLevels, m_minTemporalLevel,
                 m_spatialMeshHierarchy);
    }
    else if (m_discType == "H1H1") {
        m_disc = std::make_shared<sparseHeat::LsqSparseXtFemH1H1>
                (m_config, m_testCase,
                 m_numLevels, m_minTemporalLevel,
                 m_spatialMeshHierarchy);
    }
    else {
        std::cerr << "Unknown discretisation type!\n";
        abort();
    }
}

void sparseHeat::Solver
:: run()
{
    initialize ();
    assembleSystem();
    assembleRhs();
    solve();
}

std::pair<Vector, int> sparseHeat::Solver
:: runAndMeasurePerformanceMetrics()
{
    Vector elapsedTime(4);

    auto start = std::chrono::high_resolution_clock::now();
    initialize ();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast
            <std::chrono::milliseconds>(end - start);
    elapsedTime(0)
            = (static_cast<double>(duration.count()))/1000;

    start = std::chrono::high_resolution_clock::now();
    assembleSystem();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast
                <std::chrono::milliseconds>(end - start);
    elapsedTime(1)
                = (static_cast<double>(duration.count()))/1000;

    start = std::chrono::high_resolution_clock::now();
    assembleRhs();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast
                <std::chrono::milliseconds>(end - start);
    elapsedTime(2)
                = (static_cast<double>(duration.count()))/1000;

    int memoryUsage;
    std::tie(elapsedTime(3), memoryUsage)
            = solveAndMeasurePerformanceMetrics();

    return {elapsedTime, memoryUsage};
}

void sparseHeat::Solver
:: initialize ()
{
    // spatial FE hierarchies
    auto spatialNestedFEHierarchyForTemperature
            = m_disc->getSpatialNestedFEHierarchyForTemperature();
    auto spatialNestedFEHierarchyForHeatFlux
            = m_disc->getSpatialNestedFEHierarchyForHeatFlux();

    // solution data handler
    m_solutionHandler
            = std::make_shared<sparseHeat::SolutionHandler>
            (m_minTemporalLevel, m_maxTemporalLevel,
             spatialNestedFEHierarchyForTemperature,
             spatialNestedFEHierarchyForHeatFlux);

    // variable to store rhs data
    auto dataSize = m_solutionHandler->getDataSize();
    auto dataOffsets = evalBlockOffsets(dataSize);
    m_rhs = std::make_shared<BlockVector>(dataOffsets);
}

void sparseHeat::Solver
:: assembleSystem()
{
    m_disc->assembleSystemSubMatrices();

    if (m_linearSolver == "pardiso")
    {
        m_disc->buildSystemMatrix();
        m_systemMat = m_disc->getSystemMatrix();
    }
    else if (m_linearSolver == "cg") {
        m_disc->buildSystemMatrix();
        m_systemMat = m_disc->getSystemMatrix();
    }
}

void sparseHeat::Solver
:: assembleRhs()
{
    (*m_rhs) = 0.;
    m_disc->assembleRhs(m_rhs);
}

void sparseHeat::Solver
:: solve()
{
    // wrap rhs BlockVector as vector
    Vector wrapRhs(m_rhs->GetData(), m_rhs->Size());

    // wrap solution BlockVector as vector
    auto U = m_solutionHandler->getData();
    Vector wrapU(U->GetData(), U->Size());

    solve(wrapRhs, wrapU);
}

std::pair<double, int> sparseHeat::Solver
:: solveAndMeasurePerformanceMetrics()
{
    // wrap rhs BlockVector as vector
    Vector wrapRhs(m_rhs->GetData(), m_rhs->Size());

    // wrap solution BlockVector as vector
    auto U = m_solutionHandler->getData();
    Vector wrapU(U->GetData(), U->Size());

    return solve(wrapRhs, wrapU);
}

std::pair<double, int> sparseHeat::Solver
:: solve(const Vector& rhs, Vector& u)
{
    int memoryUsage = 0;

    auto start = std::chrono::high_resolution_clock::now();
    if (m_linearSolver == "pardiso")
    {
        setPardisoSolver();
        m_systemMat->SortColumnIndices();
        m_pardisoSolver->initialize(m_systemMat->Size(),
                                    m_systemMat->GetI(),
                                    m_systemMat->GetJ(),
                                    m_systemMat->GetData());
        m_pardisoSolver->factorize();
        m_pardisoSolver->solve(rhs.GetData(), u.GetData());
        memoryUsage = m_pardisoSolver->getMemoryUsage();
        finalizePardisoSolver();
    }
    else if (m_linearSolver == "cg")
    {
        int verbose, maxIters;
        double absTol, relTol;

        READ_CONFIG_PARAM_OR_SET_TO_DEFAULT
                (m_config, "cg_verbose", verbose, 0);
        READ_CONFIG_PARAM_OR_SET_TO_DEFAULT
                (m_config, "cg_max_iterations", maxIters, 1000);
        READ_CONFIG_PARAM_OR_SET_TO_DEFAULT
                (m_config, "cg_absolute_tolerance", absTol, 0);
        READ_CONFIG_PARAM_OR_SET_TO_DEFAULT
                (m_config, "cg_relative_tolerance", relTol, 1E-8);

        auto M = new GSSmoother(*m_systemMat);
        u = 0.;
        PCG(*m_systemMat, *M, rhs, u, verbose, maxIters, relTol, absTol);
        delete M;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast
            <std::chrono::milliseconds>(end - start);
    double elapsedTime
            = (static_cast<double>(duration.count()))/1000;

    return {elapsedTime, memoryUsage};
}

void sparseHeat::Solver
:: setPardisoSolver()
{
    int mtype = 1; // real structurally symmetric
//        int mtype = 2; // real symmetric positive definite
//        int mtype = -2; // real symmetric indefinite
    m_pardisoSolver
            = std::make_unique<PardisoSolver>(mtype);
}

void sparseHeat::Solver
:: finalizePardisoSolver() {
    m_pardisoSolver->finalize();
}

double sparseHeat::Solver
:: getMeshwidthOfFinestTemporalMesh() {
    return (m_endTime/std::pow(2, m_maxTemporalLevel));
}

double sparseHeat::Solver
:: getMeshwidthOfFinestSpatialMesh()
{
    double hMax, ddum;
    auto meshes = m_spatialMeshHierarchy->getMeshes();
    meshes[m_numLevels-1]
            ->GetCharacteristics(ddum, hMax, ddum, ddum);
    return hMax;
}

// End of file
