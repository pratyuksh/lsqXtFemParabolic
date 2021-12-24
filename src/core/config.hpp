#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <nlohmann/json.hpp>

#define READ_CONFIG_PARAM(config, key, param) \
    {                                         \
        if (config.contains(key)) {           \
            param = config[key];              \
        }                                     \
        else {                                \
            std::cout << "Key not found: "    \
                      << key << std::endl;    \
            abort();                          \
        }                                     \
    }

#define READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config,         \
                                            key,            \
                                            param,          \
                                            defaultValue)   \
    {                                                       \
        if (config.contains(key)) {                         \
            param = config[key];                            \
        }                                                   \
        else {                                              \
            param = defaultValue;                           \
        }                                                   \
    }


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
