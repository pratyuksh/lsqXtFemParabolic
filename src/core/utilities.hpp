#ifndef CORE_UTILITIES_HPP
#define CORE_UTILITIES_HPP

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <Eigen/Dense>

#include "../core/config.hpp"

//! Writes convergence behaviour of solution error to json file
void writeErrorConvergenceJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 Eigen::VectorXd &dataX1, Eigen::VectorXi &dataX2,
 Eigen::MatrixXd &dataY);

//! Writes time measurements and memory to json file
void writeMetricsJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 Eigen::VectorXd &hMax, Eigen::VectorXi &ndofs,
 Eigen::VectorXd &measuredTime, Eigen::VectorXi &memUsage);

/*
//! Writes convergence behaviour of Quantity of Interest to json file
void writeQoiConvergenceJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 Eigen::VectorXd &dataX1, Eigen::VectorXi &dataX2,
 Eigen::VectorXd &dataY1, Eigen::VectorXd &dataY2);

//! Writes time measurements to json file
void writeMeasuredTimeJsonFile
(std::string systemType, std::string solverType,
 const nlohmann::json& config,
 int lx, int lt,
 Eigen::VectorXd &measuredTime);

void write_mcmc_json_file (std::string systemType,
                           std::string solverType,
                           const nlohmann::json& config,
                           int lx, int lt,
                           Eigen::VectorXd &EQoi,
                           bool bool_mpi=false);

void write_mlmcmc_json_file (std::string systemType,
                             std::string solverType,
                             const nlohmann::json& config,
                             int lx, int lt,
                             Eigen::VectorXd &EQoi,
                             bool bool_mpi=false);
*/

#endif // CORE_UTILITIES_HPP
