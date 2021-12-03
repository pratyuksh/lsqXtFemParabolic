#include "test_cases_factory.hpp"
#include <fmt/format.h>


// Makes different test cases
std::shared_ptr<heat::TestCases>
heat::makeTestCase(const nlohmann::json& config)
{
    const std::string problem_type
            = config["problem_type"];
#ifndef NDEBUG
    std::cout << "\tFor Heat equation, run case: "
              << problem_type << std::endl;
#endif
    if (problem_type == "dummy") {
        return std::make_shared
                <heat::TestCase<Dummy>> (config);
    }
    else if (problem_type == "unitSquare_test1") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest1>> (config);
    }
    else if (problem_type == "unitSquare_test2") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest2>> (config);
    }
    else if (problem_type == "unitSquare_test3") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest3>> (config);
    }
    else if (problem_type == "unitSquare_test4") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest4>> (config);
    }
    else if (problem_type == "periodic_unitSquare_test1") {
        return std::make_shared
                <heat::TestCase<PeriodicUnitSquareTest1>> (config);
    }
    else if (problem_type == "lShaped_test1") {
        return std::make_shared
                <heat::TestCase<LShapedTest1>> (config);
    }
    else if (problem_type == "lShaped_test2") {
        return std::make_shared
                <heat::TestCase<LShapedTest2>> (config);
    }
    else if (problem_type == "lShaped_test3") {
        return std::make_shared
                <heat::TestCase<LShapedTest3>> (config);
    }

    throw std::runtime_error(fmt::format(
        "Unknown problem type for Heat equation. "
        "[{}]", problem_type));
}
