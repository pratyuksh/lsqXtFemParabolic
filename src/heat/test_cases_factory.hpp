#ifndef HEAT_TEST_CASES_FACTORY_HPP
#define HEAT_TEST_CASES_FACTORY_HPP

#include "test_cases.hpp"


namespace heat {

//! Factory of test cases for the heat equation
std::shared_ptr<heat::TestCases>
makeTestCase(const nlohmann::json& config);

}

#endif // HEAT_TEST_CASES_FACTORY_HPP
