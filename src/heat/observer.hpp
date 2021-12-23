#ifndef HEAT_OBSERVER_HPP
#define HEAT_OBSERVER_HPP

#include "mfem.hpp"

#include "../core/config.hpp"

#include "../mymfem/base_observer.hpp"
#include "../mymfem/utilities.hpp"

#include "test_cases.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "solution_handler.hpp"
//#include "functionals.hpp"


namespace heat {

/**
 * @brief Observer for the Heat equation
 *
 * Provides functions to visualize solutions, compute error
 * and write output to files.
 */
class Observer : public mymfem::BaseObserver
{
public:
    //! Default Constructor
    Observer () : BaseObserver () {}

    //! Constructor with config JSO
    //! and spatial mesh discretisation level as arguments
    Observer (const nlohmann::json&, int);

    //! Sets the test case and discretisation
    void set (std::shared_ptr<heat::TestCases>& testCase,
              std::shared_ptr<heat::LsqXtFem>& disc);

    void visualizeSolutionAtEndTime(const SolutionHandler&);

    mfem::Vector evalError(const SolutionHandler&);
    
    //! Writes the solution to output file
    void dumpSol (std::shared_ptr<mfem::GridFunction>&,
                  std::shared_ptr<mfem::GridFunction>&) const;

    /**
     * @brief Evaluates error in the natural norm
     * @param u temperature solution
     * @param q heat flux solution
     * @return error in natural norm for discrete and reference solutions
     */
    mfem::Vector evalErrorInNaturalNorm (const heat::SolutionHandler& solutionHandler) const;

    /**
     * @brief Evaluates error in the natural norm at a given time t
     * @param u temperature solution
     * @param dudt time-gradient of temperature solution
     * @param q heat flux solution
     * @param t time
     * @return error in natural norm for discrete and reference solutions
     */
    std::tuple
    <double, double, double, double, double, double>
    evalSpatialErrorOfSolutionInNaturalNorm
    (std::shared_ptr<mfem::GridFunction>& u,
     std::shared_ptr<mfem::GridFunction>& dudt,
     std::shared_ptr<mfem::GridFunction>& q,
     double t) const;

    /**
     * @brief Evaluates least-squares error
     * @param u temperature solution
     * @param q heat flux solution
     * @return least-squares error split as PDE and flux error
     */
    mfem::Vector evalErrorInLeastSquaresNorm (const heat::SolutionHandler& solutionHandler) const;

    /**
     * @brief Evaluates least-squares error at a given time t
     * @param u temperature solution
     * @param dudt time-gradient of temperature solution
     * @param q heat flux solution
     * @param t time
     * @return least-squares error split as PDE and flux error
     */
    std::tuple <double, double>
    evalSpatialErrorOfSolutionInLeastSquaresNorm
    (std::shared_ptr<mfem::GridFunction>& u,
     std::shared_ptr<mfem::GridFunction>& dudt,
     std::shared_ptr<mfem::GridFunction>& q,
     double t) const;

private:
    int m_temporalLevel;
    int m_spatialLevel;

    double m_endTime;

    bool m_boolEvalError;
    std::string m_errorType;

    std::string m_solNamePrefix;
    std::string m_solNameSuffix;

    int m_xDim;
    std::shared_ptr<heat::TestCases> m_testCase;
    std::shared_ptr<heat::LsqXtFem> m_disc;

    mutable std::unique_ptr<heat::MediumTensorCoeff> m_materialCoeff;

    mutable std::unique_ptr<heat::ExactTemperatureCoeff> m_exactTemperatureCoeff;
    mutable std::unique_ptr<heat::ExactHeatFluxCoeff> m_exactHeatFluxCoeff;

    mutable std::unique_ptr<heat::ExactTemperatureSpatialGradCoeff>
    m_exactSpatialGradientOfTemperatureCoeff;
    mutable std::unique_ptr<heat::ExactTemperatureTemporalGradCoeff>
    m_exactTemporalGradientOfTemperatureCoeff;

    mutable std::unique_ptr<heat::SourceCoeff> m_sourceCoeff;

    mfem::FiniteElementSpace* m_temporalFeSpace = nullptr;
    mfem::FiniteElementSpace* m_spatialFeSpaceForTemperature = nullptr;
    mfem::FiniteElementSpace* m_spatialFeSpaceForHeatFlux = nullptr;

    mutable std::unique_ptr<mfem::GridFunction> m_exactTemperature;
    mutable std::unique_ptr<mfem::GridFunction> m_exactHeatFlux;
};

}

#endif // HEAT_OBSERVER_HPP
