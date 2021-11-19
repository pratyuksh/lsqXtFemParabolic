#include <gtest/gtest.h>

#include "../src/core/config.hpp"
#include "../src/heat/solver.hpp"

/**
 * @brief Unit test for computing the observables for square_test1
 */
TEST(Functionals, square_test1)
{
    std::string config_file
            = "../config_files/unit_tests/heat_square_test1.json";

    auto config = getGlobalConfig(config_file);

    const int lt = config["level_t"];
    const int lx = config["level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = "../meshes/"+sub_mesh_dir;

    bool load_init_mesh = true;
    heat::Solver heat_solverFG(config, mesh_dir,
                             lx, lt, load_init_mesh);

    auto testCase = heat_solverFG.getTestCase();
    auto discr = heat_solverFG.getDiscretisation();
    auto observer = heat_solverFG.getObserver();

    auto uE_coeff
            = std::make_unique<heat::ExactTemperatureCoeff>(testCase);
    uE_coeff->SetTime(0);

    auto xVspace = discr->getXFespaces()[0];
    std::shared_ptr<GridFunction> uE
            = std::make_shared<GridFunction>(xVspace);
    uE->ProjectCoefficient(*uE_coeff);
    auto obs = observer->evalObs(uE);
    auto qoi = observer->evalQoi(uE);

    double true_obs = -8/(M_PI*M_PI);
    double true_qoi = -0.83528140191561;
    //double true_xtQoi = 1.7820427383355877;
    //std::cout << obs << "\t" << true_obs << std::endl;
    //std::cout << qoi << "\t" << true_qoi << std::endl;

    double TOL = 1E-5;
    ASSERT_LE(std::abs(obs - true_obs), TOL);
    ASSERT_LE(std::abs(qoi - true_qoi), TOL);
}

/**
 * @brief Unit test for computing the observables for periodic_square_test1
 */
TEST(Functionals, periodic_square_test1)
{
    std::string config_file
            = "../config_files/unit_tests/heat_periodic_square_test1.json";

    auto config = getGlobalConfig(config_file);

    const int lt = config["level_t"];
    const int lx = config["level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = "../meshes/"+sub_mesh_dir;

    bool load_init_mesh = true;
    heat::Solver heat_solverFG(config, mesh_dir,
                             lx, lt, load_init_mesh);

    auto testCase = heat_solverFG.getTestCase();
    auto discr = heat_solverFG.getDiscretisation();
    auto observer = heat_solverFG.getObserver();

    auto uE_coeff
            = std::make_unique<heat::ExactTemperatureCoeff>(testCase);
    uE_coeff->SetTime(1);

    auto xVspace = discr->getXFespaces()[0];
    std::shared_ptr<GridFunction> uE
            = std::make_shared<GridFunction>(xVspace);
    uE->ProjectCoefficient(*uE_coeff);
    auto obs = observer->evalObs(uE);
    auto qoi = observer->evalQoi(uE);

    double true_obs = 0;
    double true_qoi = 1.828355524808745;
    //double true_xtQoi = 1.7820427383355877;
    //std::cout << obs << "\t" << true_obs << std::endl;
    //std::cout << qoi << "\t" << true_qoi << std::endl;

    double TOL = 1E-5;
    ASSERT_LE(std::abs(obs - true_obs), TOL);
    ASSERT_LE(std::abs(qoi - true_qoi), TOL);
}
