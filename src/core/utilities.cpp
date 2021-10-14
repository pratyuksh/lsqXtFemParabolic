#include "utilities.hpp"
#include <assert.h>

namespace fs = std::filesystem;


void writeErrorConvergenceJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 Eigen::VectorXd &dataX1, Eigen::VectorXi &dataX2,
 Eigen::MatrixXd &dataY)
{
    assert(dataX1.size() == dataX2.size());

    std::string problemType = config["problem_type"];

    std::string baseOutDir = "../output";
    if (config.contains("base_out_dir")) {
        baseOutDir = config["base_out_dir"];
    }

    std::string subOutDir = "heat_"+problemType;
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

    std::string errorType = "L2H1";
    if (config.contains("error_type")) {
        errorType = config["error_type"];
    }

    std::string outfile;
    if (systemType == "heat")
    {
        outfile = outDir+"convgError"+solverType
                +"_"+discrType
                +"_"+problemType
                +"_"+errorType
                +"_deg"+std::to_string(deg)
                +".json";
    }
    std::cout << outfile << std::endl;

    // set data names
    std::string dataX1Name;
    std::string dataX2Name;
    std::vector<std::string> dataYNames;
    if (systemType == "heat")
    {
        dataX1Name = "h_max";
        dataX2Name = "ndofs";

        if (errorType == "H1") {
            dataYNames.push_back("L2");
            dataYNames.push_back("H1");
        }
        else if (errorType == "L2H1") {
            dataYNames.push_back("L2L2");
            dataYNames.push_back("L2H1");
        }
        else if (errorType == "lsq") {
            dataYNames.push_back("pde");
            dataYNames.push_back("flux");
            dataYNames.push_back("ic");
        }
        else if (errorType == "natural") {
            dataYNames.push_back("uL2H1");
            dataYNames.push_back("qL2L2");
            dataYNames.push_back("UDiv");
        }
    }

    auto json = nlohmann::json{};
    for (unsigned int i=0; i<dataX1.size(); i++) {
        json[dataX1Name][i] = dataX1(i);
        json[dataX2Name][i] = dataX2(i);
        for (unsigned int j=0; j<dataY.rows(); j++) {
            json[dataYNames[j]][i] = dataY(j,i);
        }
    }

    auto file = std::ofstream(outfile);
    assert(file.good()); // check if you can write to file.

    file << json.dump(2);
    file.close();
}

void writeMetricsJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 Eigen::VectorXd &hMax, Eigen::VectorXi &ndofs,
 Eigen::VectorXd &measuredTime, Eigen::VectorXi &memUsage)
{
    std::string problemType = config["problem_type"];

    std::string baseOutDir = "../output";
    if (config.contains("base_out_dir")) {
        baseOutDir = config["base_out_dir"];
    }

    std::string subOutDir = "heat_"+problemType;
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

    std::string outfile;
    if (systemType == "heat")
    {
        outfile = outDir+"metrics"+solverType
                +"_"+discrType
                +"_"+problemType
                +"_deg"+std::to_string(deg)
                +".json";
    }

    auto json = nlohmann::json{};
    for (unsigned int i=0; i<ndofs.size(); i++) {
        json["h_max"][i] = hMax(i);
        json["ndofs"][i] = ndofs(i);
        json["measured_time"][i] = measuredTime(i);
        json["memory_usage"][i] = memUsage(i);
    }

    auto file = std::ofstream(outfile);
    assert(file.good());
    file << json.dump(2);
    file.close();
}

/*
void writeQoiConvergenceJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 Eigen::VectorXd &dataX1, Eigen::VectorXi &dataX2,
 Eigen::VectorXd &dataY1, Eigen::VectorXd &dataY2)
{
    assert(dataX1.size() == dataX2.size());

    std::string baseOutDir = config["output_dir"];
    std::string outDir = baseOutDir+"/"+systemType;
    fs::create_directories(outDir);

    int deg = config["deg"];
    std::string baseName = config["problem_type"];

    std::string meshType = "quasiUniform";
    if (config.contains("mesh_type")) {
        meshType = config["mesh_type"];
    }

    std::string outfile;
    if (systemType == "heat")
    {
        outfile = outDir
                +"/convgQoi_"+baseName+"_"
                +meshType+"_deg"+std::to_string(deg)+"_"
                +solverType+".json";
    }

    // set data names
    std::string dataX1Name;
    std::string dataX2Name;
    std::string dataY1Name;
    std::string dataY2Name;
    if (systemType == "heat")
    {
        dataX1Name = "h_max";
        dataX2Name = "ndofs";
        dataY1Name = "qoi";
        dataY2Name = "xtQoi";
    }

    auto json = nlohmann::json{};
    for (unsigned int i=0; i<dataX1.size(); i++) {
        json[dataX1Name][i] = dataX1(i);
        json[dataX2Name][i] = dataX2(i);
        json[dataY1Name][i] = dataY1(i);
        json[dataY2Name][i] = dataY2(i);
    }

    auto file = std::ofstream(outfile);
    assert(file.good()); // check if you can write to file.

    file << json.dump(2);
    file.close();
}

void writeMeasuredtimeJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 int lx, int lt,
 Eigen::VectorXd &measuredTime)
{
    std::string baseOutDir = config["output_dir"];
    std::string outDir = baseOutDir+"/"+systemType;
    fs::create_directories(outDir);

    int deg = config["deg"];
    std::string baseName = config["problem_type"];

    std::string meshType = "quasiUniform";
    if (config.contains("mesh_type")) {
        meshType = config["mesh_type"];
    }

    auto json = nlohmann::json{};

    std::string outfile;
    if (systemType == "heat" || systemType == "mcmc_heat")
    {
        outfile = outDir
                +"/measuredTime_"+baseName+"_"
                +meshType+"_deg"+std::to_string(deg)+"_"
                +solverType+"_lx"+std::to_string(lx)
                +"_lt"+std::to_string(lt)+".json";

        json["mesh_level_x"] = lx;
        json["mesh_level_t"] = lt;
    }
    else if (systemType == "mlmcmc_heat")
    {
        int L = lx;
        int L0 = lt;

        outfile = outDir
                +"/measuredTime_"+baseName+"_"
                +meshType+"_deg"+std::to_string(deg)+"_"
                +solverType+"_L"+std::to_string(L)
                +"_L0"+std::to_string(L0)+".json";

        json["mlmcmc_max_level"] = L;
        json["mlmcmc_min_level"] = L0;
    }

    for (unsigned int i=0; i<measuredTime.size(); i++) {
        json["measured_time"][i] = measuredTime(i);
    }

    auto file = std::ofstream(outfile);
    assert(file.good()); // check if you can write to file.

    file << json.dump(2);
    file.close();
}


void write_mcmc_json_file (std::string systemType,
                           std::string solverType,
                           const nlohmann::json& config,
                           int lx, int lt,
                           Eigen::VectorXd &Eqoi,
                           bool bool_mpi)
{
    std::string base_out_dir = config["output_dir"];
    std::string out_dir = base_out_dir+"/"+systemType;
    fs::create_directories(out_dir);

    int deg = config["deg"];
    std::string base_name = config["problem_type"];

    std::string mesh_type = "quasiUniform";
    if (config.contains("mesh_type")) {
        mesh_type = config["mesh_type"];
    }

    std::string outfile;
    if (systemType == "mcmc_heat")
    {
        if (bool_mpi) {
            outfile = out_dir+"/mpi_mcmc_";
        }
        else {
            outfile = out_dir+"/mcmc_";
        }
        outfile += base_name+"_"
                +mesh_type+"_deg"+std::to_string(deg)+"_"
                +solverType+"_lx"+std::to_string(lx)
                +"_lt"+std::to_string(lt)+".json";
    }

    auto json = nlohmann::json{};
    json["mesh_level_x"] = lx;
    json["mesh_level_t"] = lt;
    for (unsigned int i=0; i<Eqoi.size(); i++) {
        json["expectation_qoi"][i] = Eqoi(i);
    }

    auto file = std::ofstream(outfile);
    assert(file.good()); // check if you can write to file.

    file << json.dump(2);
    file.close();
}

void write_mlmcmc_json_file (std::string systemType,
                             std::string solverType,
                             const nlohmann::json& config,
                             int L, int L0,
                             Eigen::VectorXd &Eqoi,
                             bool bool_mpi)
{
    std::string base_out_dir = config["output_dir"];
    std::string out_dir = base_out_dir+"/"+systemType;
    fs::create_directories(out_dir);

    int deg = config["deg"];
    std::string base_name = config["problem_type"];

    std::string mesh_type = "quasiUniform";
    if (config.contains("mesh_type")) {
        mesh_type = config["mesh_type"];
    }

    std::string outfile;
    if (systemType == "mlmcmc_heat")
    {
        if (bool_mpi) {
            outfile = out_dir+"/mpi_mlmcmc_";
        }
        else {
            outfile = out_dir+"/mlmcmc_";
        }
        outfile += base_name+"_"
                +mesh_type+"_deg"+std::to_string(deg)+"_"
                +solverType+"_L"+std::to_string(L)
                +"_L0"+std::to_string(L0)+".json";
    }

    auto json = nlohmann::json{};
    json["mlmcmc_max_level"] = L;
    json["mlmcmc_min_level"] = L0;
    for (unsigned int i=0; i<Eqoi.size(); i++) {
        json["expectation_qoi"][i] = Eqoi(i);
    }

    auto file = std::ofstream(outfile);
    assert(file.good()); // check if you can write to file.

    file << json.dump(2);
    file.close();
}*/

// End of file
