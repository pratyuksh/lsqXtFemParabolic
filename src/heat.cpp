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
    std::string problemType;
    std::string baseOutDir, subOutDir;

    int deg;
    std::string discrType, errorType;

    READ_CONFIG_PARAM(config, "problem_type", problemType);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config,
                                        "base_out_dir",
                                        baseOutDir,
                                        "../output");
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config,
                                        "sub_out_dir",
                                        subOutDir,
                                        "heat/"+problemType);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", deg, 1);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "discretisation_type",
                                        discrType, "H1Hdiv");
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "error_type",
                                        errorType, "natural");

    std::string outDir = baseOutDir+"/"+subOutDir+"/";
    fs::create_directories(outDir);

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
 Array<Vector> &elapsedTime,
 Array<int> &memoryUsage)
{
    assert(numDofs.Size() == meshSizes.Size());
    assert(elapsedTime[0].Size() == 4);

    std::string problemType;
    std::string baseOutDir, subOutDir;

    int deg;
    std::string discrType, errorType;

    READ_CONFIG_PARAM(config, "problem_type", problemType);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config,
                                        "base_out_dir",
                                        baseOutDir,
                                        "../output");
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config,
                                        "sub_out_dir",
                                        subOutDir,
                                        "heat/"+problemType);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", deg, 1);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "discretisation_type",
                                        discrType, "H1Hdiv");
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "error_type",
                                        errorType, "natural");

    std::string outDir = baseOutDir+"/"+subOutDir+"/";
    fs::create_directories(outDir);

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

        json["elapsed_time_for_initialization"][i] = elapsedTime[i].Elem(0);
        json["elapsed_time_for_system_assembly"][i] = elapsedTime[i].Elem(1);
        json["elapsed_time_for_rhs_assembly"][i] = elapsedTime[i].Elem(2);
        json["elapsed_time_for_linear_solve"][i] = elapsedTime[i].Elem(3);

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

std::tuple<int, double, double, Vector, int>
runSolverAndMeasurePerformanceMetrics(heat::Solver& solver)
{
    Vector elapsedTime;
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
    int temporalLevel, spatialLevel;
    std::string subMeshDir;

    READ_CONFIG_PARAM(config, "temporal_level", temporalLevel);
    READ_CONFIG_PARAM(config, "spatial_level", spatialLevel);
    READ_CONFIG_PARAM(config, "mesh_dir", subMeshDir);
    std::string meshDir = baseMeshDir+subMeshDir;

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
    int minTemporalLevel, minSpatialLevel, maxSpatialLevel;
    std::string subMeshDir;

    READ_CONFIG_PARAM(config, "min_temporal_level", minTemporalLevel);
    READ_CONFIG_PARAM(config, "min_spatial_level", minSpatialLevel);
    READ_CONFIG_PARAM(config, "max_spatial_level", maxSpatialLevel);
    READ_CONFIG_PARAM(config, "mesh_dir", subMeshDir);
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
    int minTemporalLevel, minSpatialLevel, maxSpatialLevel;
    std::string subMeshDir;
    int numReps;

    READ_CONFIG_PARAM(config, "min_temporal_level", minTemporalLevel);
    READ_CONFIG_PARAM(config, "min_spatial_level", minSpatialLevel);
    READ_CONFIG_PARAM(config, "max_spatial_level", maxSpatialLevel);
    READ_CONFIG_PARAM(config, "mesh_dir", subMeshDir);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT
            (config, "num_repetitions", numReps, 6);
    std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);

    int numLevels = maxSpatialLevel-minSpatialLevel+1;
    Array<double> hMax(numLevels);
    Array<int> numDofs(numLevels);
    Array<Vector> elapsedTime(numLevels);
    Array<int> memoryUsage(numLevels);

    double htMax=1, hxMax=1;
    for (int k=0; k<numLevels; k++)
    {
        int temporalLevel = minTemporalLevel+k;
        int spatialLevel = minSpatialLevel+k;

        Vector localElapsedTime;
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
            if (i == 0) {
                elapsedTime[k].SetSize(localElapsedTime.Size());
                elapsedTime[k] = 0.;
            }
            else {
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
        std::cout << "Elapsed time: ";
        elapsedTime[k].Print();
        std::cout << "Memory usage: " << memoryUsage[k] << std::endl;

        hMax[k] = (hxMax >= htMax ? hxMax : htMax);
    }
    std::cout << "Elapsed time:\n";
    for (int i=0; i<numLevels; i++) {
        elapsedTime[i].Print();
    }
    std::cout << "\nMemory usage:\n";
    memoryUsage.Print();
    std::cout << "\n#Dofs:\n";
    numDofs.Print();

    heat::writePerformanceMetricsDataToJsonFile
            (config, numDofs, hMax, elapsedTime, memoryUsage);
}


int main(int argc, char *argv[])
{   
    auto config = getGlobalConfig(argc, argv);

    std::string host, run;
    bool loadInitMesh;

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "host", host, "local");
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "run", run, "simulation");
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT
            (config, "load_init_mesh", loadInitMesh, false);

    std::string baseMeshDir;
    if (loadInitMesh) {
        baseMeshDir.assign("../meshes/");
    }
    else {
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
