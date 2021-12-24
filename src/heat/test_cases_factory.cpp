#include "test_cases_factory.hpp"
#include <fmt/format.h>


std::shared_ptr<heat::TestCases>
heat::makeTestCase(const nlohmann::json& config)
{
    std::string problemType;
    READ_CONFIG_PARAM(config, "problem_type", problemType);

#ifndef NDEBUG
    std::cout << "\tFor Heat equation, run case: "
              << problemType << std::endl;
#endif
    if (problemType == "dummy") {
        return std::make_shared
                <heat::TestCase<Dummy>> (config);
    }
    else if (problemType == "unitSquare_test1") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest1>> (config);
    }
    else if (problemType == "unitSquare_test2") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest2>> (config);
    }
    else if (problemType == "unitSquare_test3") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest3>> (config);
    }
    else if (problemType == "unitSquare_test4") {
        return std::make_shared
                <heat::TestCase<UnitSquareTest4>> (config);
    }
    else if (problemType == "periodic_unitSquare_test1") {
        return std::make_shared
                <heat::TestCase<PeriodicUnitSquareTest1>> (config);
    }
    else if (problemType == "lShaped_test1") {
        return std::make_shared
                <heat::TestCase<LShapedTest1>> (config);
    }
    else if (problemType == "lShaped_test2") {
        return std::make_shared
                <heat::TestCase<LShapedTest2>> (config);
    }
    else if (problemType == "lShaped_test3") {
        return std::make_shared
                <heat::TestCase<LShapedTest3>> (config);
    }
    else if (problemType == "unitCube_test1") {
        return std::make_shared
                <heat::TestCase<UnitCubeTest1>> (config);
    }

    throw std::runtime_error(fmt::format(
        "Unknown problem type for Heat equation. "
        "[{}]", problemType));
}
