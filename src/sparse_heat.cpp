#include "includes.hpp"
#include "core/config.hpp"
#include "core/utilities.hpp"
#include "sparse_heat/solver.hpp"

#include <iostream>
#include <Eigen/Core>


void runSolver(const nlohmann::json& config)
{
}


void runOneSimulation (const nlohmann::json config,
                       std::string baseMeshDir, 
                       bool loadInitMesh=false)
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

    if (run == "simulation") {
        runOneSimulation(config, baseMeshDir, loadInitMesh);
    }
    else if (run == "convergence") {
        // runMultipleSimulationsToTestConvergence(config, baseMeshDir, loadInitMesh);
    }

    return 0;
}


// End of file
