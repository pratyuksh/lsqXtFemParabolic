#ifndef HEAT_OBSERVER_HPP
#define HEAT_OBSERVER_HPP

#include "mfem.hpp"

#include "../core/config.hpp"

#include "../mymfem/base_observer.hpp"
#include "../mymfem/utilities.hpp"

#include "test_cases.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "functionals.hpp"


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
              std::shared_ptr<heat::LsqXtFEM>& discr);
    
    //! Writes the solution to output file
    void dumpSol (std::shared_ptr<mfem::GridFunction>&,
                  std::shared_ptr<mfem::GridFunction>&) const;
    
    /**
     * @brief Evaluates L2 and H1 error of temperature at a given time t
     * @param u temperature solution
     * @param t
     * @return L2-error and H1-error of u and reference solution
     */
    std::tuple <double, double, double, double>
    evalXH1Error (std::shared_ptr<mfem::GridFunction>& u,
                   double t) const;

    /**
     * @brief Evaluates L2L2 and L2H1 error of temperature
     * @param u Temperature solution coefficients
     * @return L2L2-error and L2H1-error of u and reference solution
     */
    std::tuple <double, double, double, double>
    evalXtL2H1Error (mfem::Vector& u) const;

    /**
     * @brief Evaluates least-squares error at a given time t
     * @param u temperature solution
     * @param dudt time-gradient of temperature solution
     * @param q heat flux solution
     * @param t time
     * @return least-squares error split as PDE and flux error
     */
    std::tuple <double, double>
    evalXLsqError (std::shared_ptr<mfem::GridFunction>& u,
                   std::shared_ptr<mfem::GridFunction>& dudt,
                   std::shared_ptr<mfem::GridFunction>& q,
                   double t) const;

    /**
     * @brief Evaluates least-squares error
     * @param u temperature solution
     * @param q heat flux solution
     * @return least-squares error split as PDE and flux error
     */
    std::tuple <double, double, double>
    evalXtLsqError (mfem::Vector& u, mfem::Vector& q) const;

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
    evalXError (std::shared_ptr<mfem::GridFunction>& u,
                std::shared_ptr<mfem::GridFunction>& dudt,
                std::shared_ptr<mfem::GridFunction>& q,
                double t) const;

    /**
     * @brief Evaluates error in the natural norm
     * @param u temperature solution
     * @param q heat flux solution
     * @return error in natural norm for discrete and reference solutions
     */
    std::tuple <double, double, double>
    evalXtError (mfem::Vector&, mfem::Vector&) const;

    double evalObs(mfem::Vector&) const;
    double evalObs(std::shared_ptr <mfem::GridFunction>&) const;

    double evalQoi(mfem::Vector&) const;
    double evalQoi(std::shared_ptr <mfem::GridFunction>&) const;

    double evalXtQoi(mfem::Vector&) const;
    
private:
    bool m_boolError;
    std::string m_solNamePrefix;
    std::string m_solNameSuffix;

    int m_xndim;
    std::shared_ptr<heat::TestCases> m_testCase;
    std::shared_ptr<heat::LsqXtFEM> m_discr;

    mutable std::unique_ptr
    <heat::MediumTensorCoeff> m_medCoeff;
    mutable std::unique_ptr
    <heat::ExactTemperatureCoeff> m_uECoeff;
    mutable std::unique_ptr<heat::ExactFluxCoeff> m_qECoeff;
    mutable std::unique_ptr
    <heat::ExactTemperatureTimeGradCoeff> m_dudtECoeff;
    mutable std::unique_ptr<heat::SourceCoeff> m_sourceCoeff;
    //mutable std::unique_ptr
    //<heat::LaplacianCoeff> m_laplacianCoeff;

    mfem::FiniteElementSpace* m_tVspace = nullptr;
    mfem::FiniteElementSpace* m_xV1space = nullptr;
    mfem::FiniteElementSpace* m_xV2space = nullptr;

    mutable std::unique_ptr<mfem::GridFunction> m_uE;
    mutable std::unique_ptr<mfem::GridFunction> m_qE;

    mutable std::unique_ptr<QoIFunctional> m_obsFn;
    mutable std::unique_ptr<QoIFunctional> m_qoiFn;
    mutable std::unique_ptr<QoIXtFunctional> m_qoiXtFn;
};

}

#endif // HEAT_OBSERVER_HPP
