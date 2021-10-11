#include "includes.hpp"
#include "core/config.hpp"
#include "core/utilities.hpp"
#include "heat/solver.hpp"

#include <iostream>
#include <Eigen/Core>
#include <chrono>


void runFG (const nlohmann::json& config,
            std::string baseMeshDir, bool loadInitMesh=false)
{
    const int lt = config["level_t"];
    const int lx = config["level_x"];
    std::string subMeshDir = config["mesh_dir"];
    const std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);
    heat::Solver heatSolverFG(config, testCase,
                              meshDir, lx, lt, loadInitMesh);

    int ndofs;
    double htMax, hxMax;
    Eigen::VectorXd errSol;
    std::tie (ndofs, htMax, hxMax, errSol) = heatSolverFG();
    std::cout << "\n\nError: "
              << errSol.transpose() << std::endl;
    std::cout << "tMesh size: " << htMax << std::endl;
    std::cout << "xMesh size: " << hxMax << std::endl;
    std::cout << "#Dofs: " << ndofs << std::endl;
}

void convergenceFG (const nlohmann::json& config,
                    std::string baseMeshDir, bool loadInitMesh=false)
{
    const int Lt0 = config["min_level_t"];
    const int Lx0 = config["min_level_x"];
    const int Lx = config["max_level_x"];
    std::string subMeshDir = config["mesh_dir"];
    const std::string meshDir = baseMeshDir+subMeshDir;

    int num_levels = Lx-Lx0+1;
    Eigen::VectorXd h_max(num_levels);
    Eigen::VectorXi ndofs(num_levels);
    Eigen::MatrixXd errSol(2,num_levels);
    errSol.setZero();

    heat::Solver *heatSolverFG = nullptr;
    double htMax, hxMax;
    Eigen::VectorXd errSolBuf;
    for (int k=0; k<num_levels; k++)
    {
        int lt = Lt0+k;
        int lx = Lx0+k;
        heatSolverFG = new heat::Solver(config, meshDir,
                                        lx, lt, loadInitMesh);
        std::tie (ndofs(k), htMax,
                  hxMax, errSolBuf) = (*heatSolverFG)();
        if(k == 0) {
            errSol.resize(errSolBuf.size(), num_levels);
        }
        errSol.col(k) = errSolBuf;
        std::cout << "\nLevels: " << lt << ", "
                  << lx << std::endl;
        std::cout << "tMesh size: " << htMax << std::endl;
        std::cout << "xMesh size: " << hxMax << std::endl;
        std::cout << "#Dofs: " << ndofs(k) << std::endl;
        std::cout << "Error: "
                  << errSolBuf.transpose() << std::endl;
        delete heatSolverFG;

        h_max(k) = (hxMax >= htMax ? hxMax : htMax);
    }
    std::cout << "\n\nError:\n"
              << errSol.transpose() << std::endl;
    std::cout << "\n#Dofs: "
              << ndofs.transpose() << std::endl;

    // write convergence results to json file
    writeErrorConvergenceJsonFile ("heat", "FG",
                                   config, h_max, ndofs, errSol);
}

void measureMetricsFG (const nlohmann::json& config,
                       std::string baseMeshDir, bool loadInitMesh=false)
{
    const int Lt0 = config["min_level_t"];
    const int Lx0 = config["min_level_x"];
    const int Lx = config["max_level_x"];
    std::string subMeshDir = config["mesh_dir"];
    const std::string meshDir = baseMeshDir+subMeshDir;

    int numReps = 1;
    if (config.contains("num_reps")) {
        numReps = config["num_reps"];
    }

    int num_levels = Lx-Lx0+1;
    Eigen::VectorXd h_max(num_levels);
    Eigen::VectorXi ndofs(num_levels);
    Eigen::VectorXd elapsedTime(num_levels);
    Eigen::VectorXi memUsage(num_levels);

    heat::Solver *solver = nullptr;
    double htMax, hxMax;
    for (int k=0; k<num_levels; k++)
    {
        int lt = Lt0+k;
        int lx = Lx0+k;
        solver = new heat::Solver(config, meshDir, lx, lt, loadInitMesh);

        std::tie (ndofs(k), htMax, hxMax,
                  elapsedTime(k), memUsage(k)) = solver->measure(numReps);
        std::cout << "\nLevels: " << lt << ", "
                  << lx << std::endl;
        std::cout << "tMesh size: " << htMax << std::endl;
        std::cout << "xMesh size: " << hxMax << std::endl;
        std::cout << "#Dofs: " << ndofs(k) << std::endl;
        std::cout << "Solve time: " << elapsedTime(k) << " seconds" << std::endl;
        std::cout << "Memory usage: " << memUsage(k) << " kB\n" << std::endl;
        delete solver;

        h_max(k) = (hxMax >= htMax ? hxMax : htMax);
    }
    std::cout << "\n\n#Dofs: "
              << ndofs.transpose() << std::endl;
    std::cout << "\nSolve time (s):\n"
              << elapsedTime.transpose() << std::endl;
    std::cout << "\nMemory usage (kB):\n"
              << memUsage.transpose() << std::endl;


    // write time measurement results to json file
    writeMetricsJsonFile("heat", "FG",
                         config, h_max, ndofs, elapsedTime, memUsage);
}


int main(int argc, char *argv[])
{   
    // Read config json
    auto config = getGlobalConfig(argc, argv);
    const std::string host = config["host"];
    const std::string run = config["run"];
    std::string baseMeshDir;

    // check if an initial mesh needs to be loaded
    bool loadInitMesh = false;
    if (config.contains("load_init_mesh")) {
        loadInitMesh = config["load_init_mesh"];
    }

    if (loadInitMesh) {
        baseMeshDir.assign("../meshes/");
    } else {
        if (host == "local") {
            baseMeshDir.assign(localBaseMeshDir);
        }
        else if (host == "cluster") {
            baseMeshDir.assign(clusterBaseMeshDir);
        }
    }

    if (run == "runFG") {
        runFG(config, baseMeshDir, loadInitMesh);
    }
    else if (run == "convergenceFG") {
        convergenceFG(config, baseMeshDir, loadInitMesh);
    }
    else if (run == "measureFG") {
        measureMetricsFG(config, baseMeshDir, loadInitMesh);
    }

    return 0;
}


// End of file
