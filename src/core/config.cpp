#include "config.hpp"
#include <iostream>
#include <fstream>


//! Reads the second argument from the command line as config file name,
//! else throws an error.
//! Creates the JSON config object from the config file name
nlohmann::json getGlobalConfig(int argc, char* const argv[])
{
    nlohmann::json config;
    if (argc == 2) {
        std::string fileName = argv[1];
        config = getGlobalConfig (fileName);
	}
	else {
		std::cout << "\nError: Wrong number of arguments.\n";
	}
	
	return config;
}

//! Creates the JSON config object from the config file name
nlohmann::json getGlobalConfig(const std::string& fileName)
{
    std::ifstream file(fileName);
    nlohmann::json config;
    file >> config;
    return config;
}
