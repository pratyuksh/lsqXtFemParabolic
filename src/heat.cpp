#include "includes.hpp"
#include "core/config.hpp"
#include "heat/solver.hpp"
#include "heat/observer.hpp"

#include <iostream>
#include <filesystem>

using namespace mfem;
namespace fs = std::filesystem;

namespace heat
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

    std::string subOutDir = "heat/"+problemType;
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

void writePerformanceMetricsDataToJsonFile
(const nlohmann::json& config,
 Array<int> &numDofs,
 Array<double> &meshSizes,
 Array<double> &elapsedTime,
 Array<int> &memoryUsage)
{
    assert(numDofs.Size() == meshSizes.Size());

    std::string problemType = config["problem_type"];

    std::string baseOutDir = "../output";
    if (config.contains("base_out_dir")) {
        baseOutDir = config["base_out_dir"];
    }

    std::string subOutDir = "heat/"+problemType;
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

    std::string outfile = outDir+"metrics"
            +"_"+discrType
            +"_"+problemType
            +"_"+errorType
            +"_deg"+std::to_string(deg)
            +".json";;
    std::cout << outfile << std::endl;

    auto json = nlohmann::json{};
    for (int i=0; i<numDofs.Size(); i++) {
        json["ndofs"][i] = numDofs[i];
        json["h_max"][i] = meshSizes[i];
        json["elapsed_time"][i] = elapsedTime[i];
        json["memory_usage"][i] = memoryUsage[i];
    }

    auto file = std::ofstream(outfile);
    assert(file.good());
    file << json.dump(2);
    file.close();
}

}


std::tuple<int, double, double, Vector>
runSolver(heat::Solver& solver,
          heat::Observer& observer)
{
    solver.run();
    auto solutionHandler = solver.getSolutionHandler();
    auto testCase = solver.getTestCase();
    auto disc = solver.getDiscretisation();

    observer.set(testCase, disc);
    observer.visualizeSolutionAtEndTime(*solutionHandler);
    auto solutionError = observer.evalError(*solutionHandler);

    double htMax, hxMax;
    std::tie(htMax, hxMax) = solver.getMeshwidths();

    return {solver.getNumDofs(), htMax, hxMax, solutionError};
}

std::tuple<int, double, double, double, int>
runSolverAndMeasurePerformanceMetrics(heat::Solver& solver)
{
    double elapsedTime;
    int memoryUsage;
    std::tie(elapsedTime, memoryUsage)
            = solver.runAndMeasurePerformanceMetrics();

    double htMax, hxMax;
    std::tie(htMax, hxMax) = solver.getMeshwidths();

    return {solver.getNumDofs(), htMax, hxMax,
                elapsedTime, memoryUsage};
}

void runOneSimulation (const nlohmann::json config,
                       std::string baseMeshDir,
                       bool loadInitMesh=false)
{
    const int temporalLevel = config["level_t"];
    const int spatialLevel = config["level_x"];
    std::string subMeshDir = config["mesh_dir"];
    const std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);

    heat::Solver solver(config,
                        testCase,
                        meshDir,
                        spatialLevel,
                        temporalLevel,
                        loadInitMesh);

    heat::Observer observer(config, spatialLevel);

    int numDofs;
    double htMax, hxMax;
    Vector solutionError;
//    std::tie (numDofs, htMax, hxMax, solutionError) = solver();
    std::tie (numDofs, htMax, hxMax, solutionError)
            = runSolver(solver, observer);

    std::cout << "\n\nSolution Error: ";
    solutionError.Print();
    std::cout << "tMesh size: " << htMax << std::endl;
    std::cout << "xMesh size: " << hxMax << std::endl;
    std::cout << "#Dofs: " << numDofs << std::endl;
}

void runMultipleSimulationsToTestConvergence (const nlohmann::json config,
                                              std::string baseMeshDir,
                                              bool loadInitMesh=false)
{
    int minTemporalLevel = config["min_level_t"];
    int minSpatialLevel = config["min_level_x"];
    int maxSpatialLevel = config["max_level_x"];
    std::string subMeshDir = config["mesh_dir"];
    std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);

    int numLevels = maxSpatialLevel-minSpatialLevel+1;
    Array<double> hMax(numLevels);
    Array<int> numDofs(numLevels);
    Array<Vector> solutionError(numLevels);

    double htMax, hxMax;
    Vector errSolBuf;
    for (int k=0; k<numLevels; k++)
    {
        int temporalLevel = minTemporalLevel+k;
        int spatialLevel = minSpatialLevel+k;

        heat::Solver solver (config,
                             testCase,
                             meshDir,
                             spatialLevel,
                             temporalLevel,
                             loadInitMesh);

        heat::Observer observer (config, spatialLevel);

        std::tie (numDofs[k], htMax, hxMax, solutionError[k])
                = runSolver(solver, observer);

        std::cout << "\nLevels: "
                  << temporalLevel << ", "
                  << spatialLevel << std::endl;
        std::cout << "tMesh size: " << htMax << std::endl;
        std::cout << "xMesh size: " << hxMax << std::endl;
        std::cout << "#Dofs: " << numDofs[k] << std::endl;
        std::cout << "Error: ";
        solutionError[k].Print();

        hMax[k] = (hxMax >= htMax ? hxMax : htMax);
    }
    std::cout << "\n\nError:\n";
    for (int i=0; i<numLevels; i++) {
        solutionError[i].Print();
    }
    std::cout << "\n#Dofs:\n";
    numDofs.Print();

    // write convergence results to json file
    heat::writeErrorConvergenceDataToJsonFile(config, numDofs,
                                              hMax, solutionError);
}

void runMultipleSimulationsToMeasurePerformanceMetrics
(const nlohmann::json config,
 std::string baseMeshDir,
 bool loadInitMesh=false)
{
    int minTemporalLevel = config["min_level_t"];
    int minSpatialLevel = config["min_level_x"];
    int maxSpatialLevel = config["max_level_x"];
    std::string subMeshDir = config["mesh_dir"];
    std::string meshDir = baseMeshDir+subMeshDir;

    int numReps = 6;
    if (config.contains("num_repetitions")) {
        numReps = config["num_repetitions"];
    }

    auto testCase = heat::makeTestCase(config);

    int numLevels = maxSpatialLevel-minSpatialLevel+1;
    Array<double> hMax(numLevels);
    Array<int> numDofs(numLevels);
    Array<double> elapsedTime(numLevels);
    Array<int> memoryUsage(numLevels);

    double htMax=1, hxMax=1;
    for (int k=0; k<numLevels; k++)
    {
        int temporalLevel = minTemporalLevel+k;
        int spatialLevel = minSpatialLevel+k;

        double localElapsedTime;
        elapsedTime[k] = 0;
        for (int i=0; i<numReps; i++)
        {
            heat::Solver solver (config,
                                 testCase,
                                 meshDir,
                                 spatialLevel,
                                 temporalLevel,
                                 loadInitMesh);

            std::tie (numDofs[k], htMax, hxMax,
                      localElapsedTime, memoryUsage[k])
                    = runSolverAndMeasurePerformanceMetrics(solver);
            if (i > 0) {
                elapsedTime[k] += localElapsedTime;
            }
        }
        elapsedTime[k] /= (numReps-1);

        std::cout << "\nLevels: "
                  << temporalLevel << ", "
                  << spatialLevel << std::endl;
        std::cout << "tMesh size: " << htMax << std::endl;
        std::cout << "xMesh size: " << hxMax << std::endl;
        std::cout << "#Dofs: " << numDofs[k] << std::endl;
        std::cout << "Elapsed time: " << elapsedTime[k] << std::endl;
        std::cout << "Memory usage: " << memoryUsage[k] << std::endl;

        hMax[k] = (hxMax >= htMax ? hxMax : htMax);
    }
    std::cout << "\n\nElapsed time:\n";
    elapsedTime.Print();
    std::cout << "\nMemory usage:\n";
    memoryUsage.Print();
    std::cout << "\n#Dofs:\n";
    numDofs.Print();

    heat::writePerformanceMetricsDataToJsonFile
            (config, numDofs, hMax, elapsedTime, memoryUsage);
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
    else if (run == "measure") {
        runMultipleSimulationsToMeasurePerformanceMetrics
                (std::move(config), std::move(baseMeshDir), loadInitMesh);
    }

    return 0;
}


// End of file
