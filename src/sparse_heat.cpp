#include "includes.hpp"
#include "core/config.hpp"
#include "core/utilities.hpp"
#include "sparse_heat/solver.hpp"
#include "sparse_heat/observer.hpp"

#include <iostream>

using namespace mfem;


std::tuple<int, double, double, Vector>
runSolver(sparseHeat::Solver& solver,
          sparseHeat::Observer& observer)
{
    solver.run();
    auto solutionHandler = solver.getSolutionHandler();
    auto testCase = solver.getTestCase();
    auto disc = solver.getDiscretisation();

    observer.set(testCase, disc);
    auto solutionError = observer.evalError(*solutionHandler);

    double htMax = solver.getMeshwidthOfFinestTemporalMesh();
    double hxMax = solver.getMeshwidthOfFinestSpatialMesh();

    return {solver.getNumDofs(), htMax, hxMax, solutionError};
}


void runOneSimulation (const nlohmann::json config,
                       std::string baseMeshDir, 
                       bool loadInitMesh=false)
{
    int deg = config["deg"];
    assert(deg == 1);

    auto testCase = heat::makeTestCase(config);

    std::string subMeshDir = config["mesh_dir"];
    const std::string meshDir = baseMeshDir+subMeshDir;

    int numLevels = 1;
    if (config.contains("num_levels")) {
        numLevels = config["num_levels"];
    }
    int minSpatialLevel = 1;
    if (config.contains("min_spatial_level")) {
        minSpatialLevel = config["min_spatial_level"];
    }

    int minTemporalLevel = 1;
    if (config.contains("min_temporal_level")) {
        minTemporalLevel = config["min_temporal_level"];
    }

    sparseHeat::Solver solver(config,
                              testCase,
                              meshDir,
                              numLevels,
                              minSpatialLevel,
                              minTemporalLevel,
                              loadInitMesh);

    sparseHeat::Observer observer(config, numLevels, minTemporalLevel);

    int numDofs;
    double htMax, hxMax;
    Vector solutionError;
    std::tie(numDofs, htMax, hxMax, solutionError)
            = runSolver(solver, observer);

    std::cout << "\nSolution error: ";
    solutionError.Print();
    std::cout << "Temporal mesh size: " << htMax << std::endl;
    std::cout << "Spatial mesh size: " << hxMax << std::endl;
    std::cout << "#Dofs: " << numDofs << std::endl;
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
        runOneSimulation(std::move(config),
                         std::move(baseMeshDir),
                         loadInitMesh);
    }
    else if (run == "convergence") {
        // runMultipleSimulationsToTestConvergence(config, baseMeshDir, loadInitMesh);
    }

    return 0;
}


// End of file
