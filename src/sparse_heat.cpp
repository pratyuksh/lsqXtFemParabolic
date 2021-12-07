#include "includes.hpp"
#include "core/config.hpp"
#include "core/utilities.hpp"
#include "sparse_heat/solver.hpp"
#include "sparse_heat/observer.hpp"

#include <iostream>

using namespace mfem;
namespace fs = std::filesystem;

namespace sparseHeat
{

void writeErrorConvergenceDataToJsonFile
(const nlohmann::json& config,
 Array<int> &numDofs,
 Array<double> &meshSizes,
 Array<Vector> &solutionError)
{
    assert(numDofs.Size() == meshSizes.Size());

    std::string problemType = config["problem_type"];

    std::string baseOutDir = "../output";
    if (config.contains("base_out_dir")) {
        baseOutDir = config["base_out_dir"];
    }

    std::string subOutDir = "spars_heat/"+problemType;
    if (config.contains("sub_out_dir")) {
        subOutDir = config["sub_out_dir"];
    }

    std::string outDir = baseOutDir+"/"+subOutDir+"/";
    fs::create_directories(outDir);

    int deg = config["deg"];

    std::string discrType = "H1Hdiv";
    if (config.contains("discretisation_type")) {
        discrType = config["discretisation_type"];
    }

    std::string errorType = "natural";
    if (config.contains("error_type")) {
        errorType = config["error_type"];
    }

    std::string outfile = outDir+"convgError"
            +"_"+discrType
            +"_"+problemType
            +"_"+errorType
            +"_deg"+std::to_string(deg)
            +".json";;
    std::cout << outfile << std::endl;

    // set data names
    std::vector<std::string> solutionErrorNames;
    if (errorType == "natural") {
            solutionErrorNames.push_back("uL2H1");
            solutionErrorNames.push_back("qL2L2");
            solutionErrorNames.push_back("UDiv");
        }
    else if (errorType == "lsq") {
        solutionErrorNames.push_back("pde");
        solutionErrorNames.push_back("flux");
        solutionErrorNames.push_back("ic");
    }

    auto json = nlohmann::json{};
    for (int i=0; i<numDofs.Size(); i++) {
        json["ndofs"][i] = numDofs[i];
        json["h_max"][i] = meshSizes[i];
        for (int j=0; j<solutionError[i].Size(); j++) {
            json[solutionErrorNames[j]][i] = (solutionError[i])(j);
        }
    }

    auto file = std::ofstream(outfile);
    assert(file.good());
    file << json.dump(2);
    file.close();
}

}


std::tuple<int, double, double, Vector>
runSolver(sparseHeat::Solver& solver,
          sparseHeat::Observer& observer)
{
    solver.run();
    auto solutionHandler = solver.getSolutionHandler();
    auto testCase = solver.getTestCase();
    auto disc = solver.getDiscretisation();

    observer.set(testCase, disc);
    observer.visualizeSolutionAtEndTime(*solutionHandler);
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


void runMultipleSimulationsToTestConvergence (const nlohmann::json config,
                                              std::string baseMeshDir,
                                              bool loadInitMesh=false)
{
    int deg = config["deg"];
    assert(deg == 1);

    auto testCase = heat::makeTestCase(config);

    std::string subMeshDir = config["mesh_dir"];
    const std::string meshDir = baseMeshDir+subMeshDir;

    int minNumLevels = 1;
    if (config.contains("min_num_levels")) {
        minNumLevels = config["min_num_levels"];
    }

    int maxNumLevels = 1;
    if (config.contains("max_num_levels")) {
        maxNumLevels = config["max_num_levels"];
    }

    int minSpatialLevel = 1;
    if (config.contains("min_spatial_level")) {
        minSpatialLevel = config["min_spatial_level"];
    }

    int minTemporalLevel = 1;
    if (config.contains("min_temporal_level")) {
        minTemporalLevel = config["min_temporal_level"];
    }

    Array<int> numDofs(maxNumLevels-minNumLevels+1);
    Array<double> htMax(numDofs.Size());
    Array<double> hxMax(numDofs.Size());
    Array<Vector> solutionError(numDofs.Size());

    int count = 0;
    for (int numLevels = minNumLevels;
         numLevels <= maxNumLevels; numLevels++)
    {
        sparseHeat::Solver solver(config,
                                  testCase,
                                  meshDir,
                                  numLevels,
                                  minSpatialLevel,
                                  minTemporalLevel,
                                  loadInitMesh);

        sparseHeat::Observer observer(config, numLevels, minTemporalLevel);

        std::tie(numDofs[count], htMax[count], hxMax[count], solutionError[count])
                = runSolver(solver, observer);

        count++;
    }

    std::cout << "\nTemporal mesh sizes: "; htMax.Print();
    std::cout << "Spatial mesh sizes: "; hxMax.Print();
    std::cout << "#Dofs: "; numDofs.Print();
    std::cout << "Solution errors:\n";
    for (int i=0; i<maxNumLevels-minNumLevels+1; i++) {
        solutionError[i].Print();
    }

    Array<double> meshSizes(numDofs.Size());
    for (int i=0; i<maxNumLevels-minNumLevels+1; i++) {
        meshSizes[i] = hxMax[i] > htMax[i] ? hxMax[i] : htMax[i];
    }

    sparseHeat::writeErrorConvergenceDataToJsonFile(config,
                                                    numDofs,
                                                    meshSizes,
                                                    solutionError);
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
        runMultipleSimulationsToTestConvergence(std::move(config),
                                                std::move(baseMeshDir),
                                                loadInitMesh);
    }

    return 0;
}


// End of file
