#ifndef HEAT_FUNCTIONALS_HPP
#define HEAT_FUNCTIONALS_HPP

#include "mfem.hpp"

namespace heat {

/**
 * @brief Quantity of interest LF Integrator; (x^a, fluxQ)
 */
class QoILFIntegrator
        : public mfem::LinearFormIntegrator
{
public:
    //! Constructor with power exponent provided as argument
    QoILFIntegrator(double a) : m_pow(a) {}

    //! Assembles the integrator on a mesh element
    void AssembleRHSElementVect (const mfem::FiniteElement &,
                                 mfem::ElementTransformation &,
                                 mfem::Vector &) override;

private:
    double m_pow;
#ifndef MFEM_THREAD_SAFE
    DenseMatrix dshape;
#endif
};

/**
 * @brief Quantity of interest functional; (x^a, fluxQ)
 */
class QoIFunctional
{
public:
    //! Constructor with FE space passed as an argument
    QoIFunctional(mfem::FiniteElementSpace *fes)
        : m_fes(fes) {
        init();
    }

    //! Constructor with FE space and power exponent passed as argument
    QoIFunctional(mfem::FiniteElementSpace *fes, double a)
        : m_pow(a), m_fes(fes) {
        init();
    }

    //! Default destructor
    ~QoIFunctional () {}

    //! Initializes the integrator
    void init() const;

    //! Evalautes the QoI for a function passed as an MFEM GridFunction
    double operator()(const mfem::GridFunction& f) const;

    //! Evalautes the QoI for a function passed as a vector of coefficients
    double operator()(const mfem::Vector& f) const;

private:
    double m_pow = 1;
    mfem::FiniteElementSpace *m_fes = nullptr;
    mutable mfem::Vector m_bObs;
};

/**
 * @brief Quantity of interest xt-functional; (x^a, fluxQ)
 */
class QoIXtFunctional
{
public:
    //! Constructor with space-time FE spaces passed as arguments
    QoIXtFunctional(mfem::FiniteElementSpace *tFes,
                    mfem::FiniteElementSpace *xFes)
        : m_tFes(tFes), m_xFes(xFes) {
        init();
    }

    //! Constructor with space-time FE spaces
    //! and power exponent passed as arguments
    QoIXtFunctional(mfem::FiniteElementSpace *tFes,
                    mfem::FiniteElementSpace *xFes,
                    double a)
        : m_pow(a), m_tFes(tFes), m_xFes(xFes) {
        init();
    }

    //! Default destructor
    ~QoIXtFunctional () {}

    //! Initializes the integrator
    void init() const;

    //! Evaluate the QoI
    double operator()(const mfem::Vector& f) const;

private:
    void add_vector(const mfem::Array<int> &vdofs,
                    const mfem::Vector& elvec,
                    mfem::Vector& b) const;

    double m_pow = 1;
    mfem::FiniteElementSpace *m_tFes = nullptr;
    mfem::FiniteElementSpace *m_xFes = nullptr;
    mutable mfem::Vector m_bObs;
};

}

#endif // HEAT_FUNCTIONALS_HPP
