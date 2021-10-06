#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <nlohmann/json.hpp>

/**
 * @brief Generates a JSON config object.
 * @param argc number of arguments
 * @param argv argument values
 * @return JSON config object
 */
nlohmann::json getGlobalConfig(int argc, char* const argv[]);

/**
 * @brief Generates a JSON config object.
 * @param fileName name of config file
 * @return JSON config object
 */
nlohmann::json getGlobalConfig(const std::string& fileName);

#endif // CONFIG_HPP
