#include "includes.hpp"
#include "core/config.hpp"
#include "sparse_heat/solver.hpp"
#include "sparse_heat/observer.hpp"

#include <iostream>
#include <filesystem>

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
                                        "sparse_heat/"+problemType);

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
                                        "sparse_heat/"+problemType);

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

std::tuple<int, double, double, Vector, int>
runSolverAndMeasurePerformanceMetrics(sparseHeat::Solver& solver)
{
    Vector elapsedTime;
    int memoryUsage;
    std::tie(elapsedTime, memoryUsage)
            = solver.runAndMeasurePerformanceMetrics();

    double htMax = solver.getMeshwidthOfFinestTemporalMesh();
    double hxMax = solver.getMeshwidthOfFinestSpatialMesh();

    return {solver.getNumDofs(), htMax, hxMax,
                elapsedTime, memoryUsage};
}


void runOneSimulation (const nlohmann::json config,
                       std::string baseMeshDir, 
                       bool loadInitMesh=false)
{
    int deg;
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", deg, 1);
    assert(deg == 1);

    int numLevels, minSpatialLevel, minTemporalLevel;
    std::string subMeshDir;

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "num_levels",
                                        numLevels, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_spatial_level",
                                        minSpatialLevel, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_temporal_level",
                                        minTemporalLevel, 1);
    READ_CONFIG_PARAM(config, "mesh_dir", subMeshDir);
    std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);

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
    int deg;
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", deg, 1);
    assert(deg == 1);

    int minNumLevels, maxNumLevels;
    int minSpatialLevel, minTemporalLevel;
    std::string subMeshDir;

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_num_levels",
                                        minNumLevels, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "max_num_levels",
                                        maxNumLevels, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_spatial_level",
                                        minSpatialLevel, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_temporal_level",
                                        minTemporalLevel, 1);
    READ_CONFIG_PARAM(config, "mesh_dir", subMeshDir);

    std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);

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

        std::tie(numDofs[count], htMax[count], hxMax[count],
                 solutionError[count])
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

void runMultipleSimulationsToMeasurePerformanceMetrics
(const nlohmann::json config,
 std::string baseMeshDir,
 bool loadInitMesh=false)
{
    int deg;
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", deg, 1);
    assert(deg == 1);

    int minNumLevels, maxNumLevels;
    int minSpatialLevel, minTemporalLevel;
    std::string subMeshDir;
    int numReps;

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_num_levels",
                                        minNumLevels, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "max_num_levels",
                                        maxNumLevels, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_spatial_level",
                                        minSpatialLevel, 1);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "min_temporal_level",
                                        minTemporalLevel, 1);
    READ_CONFIG_PARAM(config, "mesh_dir", subMeshDir);
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "num_repetitions",
                                        numReps, 6);

    std::string meshDir = baseMeshDir+subMeshDir;

    auto testCase = heat::makeTestCase(config);

    Array<int> numDofs(maxNumLevels-minNumLevels+1);
    Array<double> htMax(numDofs.Size());
    Array<double> hxMax(numDofs.Size());
    Array<Vector> elapsedTime(numDofs.Size());
    Array<int> memoryUsage(numDofs.Size());

    int count = 0;
    for (int numLevels = minNumLevels;
         numLevels <= maxNumLevels; numLevels++)
    {
        Vector localElapsedTime;
        for (int i=0; i<numReps; i++)
        {
            sparseHeat::Solver solver(config,
                                      testCase,
                                      meshDir,
                                      numLevels,
                                      minSpatialLevel,
                                      minTemporalLevel,
                                      loadInitMesh);

            std::tie(numDofs[count], htMax[count], hxMax[count],
                     localElapsedTime, memoryUsage[count])
                    = runSolverAndMeasurePerformanceMetrics(solver);
            if (i == 0) {
                elapsedTime[count].SetSize(localElapsedTime.Size());
                elapsedTime[count] = 0.;
            }
            else {
                elapsedTime[count] += localElapsedTime;
            }
        }
        elapsedTime[count] /= (numReps-1);

        count++;
    }

    std::cout << "\nTemporal mesh sizes: "; htMax.Print();
    std::cout << "Spatial mesh sizes: "; hxMax.Print();
    std::cout << "#Dofs: "; numDofs.Print();
    std::cout << "Elapsed time:\n";
    for (int i=0; i<maxNumLevels-minNumLevels+1; i++) {
        elapsedTime[i].Print();
    }
    std::cout << "\nMemory usage:\n";
    memoryUsage.Print();

    Array<double> meshSizes(numDofs.Size());
    for (int i=0; i<maxNumLevels-minNumLevels+1; i++) {
        meshSizes[i] = hxMax[i] > htMax[i] ? hxMax[i] : htMax[i];
    }

    sparseHeat::writePerformanceMetricsDataToJsonFile(config,
                                                      numDofs,
                                                      meshSizes,
                                                      elapsedTime,
                                                      memoryUsage);
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
                (std::move(config),
                 std::move(baseMeshDir),
                 loadInitMesh);
    }

    return 0;
}


// End of file
