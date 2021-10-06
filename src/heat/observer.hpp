#ifndef HEAT_OBSERVER_HPP
#define HEAT_OBSERVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../core/config.hpp"

#include "../mymfem/base_observer.hpp"
#include "../mymfem/utilities.hpp"

#include "test_cases.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "functionals.hpp"


namespace heat {

// Observer for Heat system
class Observer : public mymfem::BaseObserver
{
public:
    Observer () : BaseObserver () {}

    Observer (const nlohmann::json&, int);

    void set (std::shared_ptr<heat::TestCases>& testCase,
              std::shared_ptr<heat::LsqXtFEM>& discr);
    
    void dumpSol (std::shared_ptr<GridFunction>&,
                   std::shared_ptr<GridFunction>&) const;
    
    //! Computes H1 error at a given time t
    std::tuple <double, double, double, double>
    evalXH1Error (std::shared_ptr<GridFunction>&,
                   double t) const;

    //! Computes L2H1 error
    std::tuple <double, double, double, double>
    evalXtL2H1Error (Vector&) const;

    //! Computes least-squares error at a given time t
    std::tuple <double, double>
    evalXLsqError (std::shared_ptr<GridFunction>&,
                    std::shared_ptr<GridFunction>&,
                    std::shared_ptr<GridFunction>&,
                    double t) const;

    //! Computes least-squares error
    std::tuple <double, double, double>
    evalXtLsqError (Vector&, Vector&) const;

    //! Computes the error in the natural norm
    //! at a given time t
    std::tuple
    <double, double, double, double, double, double>
    evalXError (std::shared_ptr<GridFunction>&,
                 std::shared_ptr<GridFunction>&,
                 std::shared_ptr<GridFunction>&,
                 double t) const;

    //! Computes the error in the natural norm
    std::tuple <double, double, double>
    evalXtError (Vector&, Vector&) const;

    double evalObs(Vector&) const;
    double evalObs(std::shared_ptr <GridFunction>&) const;

    double evalQoi(Vector&) const;
    double evalQoi(std::shared_ptr <GridFunction>&) const;

    double evalXtQoi(Vector&) const;
    
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

    FiniteElementSpace* m_tVspace = nullptr;
    FiniteElementSpace* m_xV1space = nullptr;
    FiniteElementSpace* m_xV2space = nullptr;

    mutable std::unique_ptr<GridFunction> m_uE;
    mutable std::unique_ptr<GridFunction> m_qE;

    mutable std::unique_ptr<QoIFunctional> m_obsFn;
    mutable std::unique_ptr<QoIFunctional> m_qoiFn;
    mutable std::unique_ptr<QoIXtFunctional> m_qoiXtFn;
};

}

#endif // HEAT_OBSERVER_HPP
