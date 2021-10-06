#ifndef HEAT_FUNCTIONALS_HPP
#define HEAT_FUNCTIONALS_HPP

#include "mfem.hpp"
using namespace mfem;

namespace heat {

// Quantity of interest LF Integrator; (x^a, fluxQ)
class QoILFIntegrator
        : public LinearFormIntegrator
{
public:
    QoILFIntegrator(double a) : m_pow(a) {}

    void AssembleRHSElementVect (const FiniteElement &,
                                 ElementTransformation &,
                                 Vector &) override;

private:
    double m_pow;
#ifndef MFEM_THREAD_SAFE
    DenseMatrix dshape;
#endif
};

// Quantity of interest functional; (x^a, fluxQ)
class QoIFunctional
{
public:

    QoIFunctional(FiniteElementSpace *fes)
        : m_fes(fes) {
        init();
    }

    QoIFunctional(FiniteElementSpace *fes, double a)
        : m_pow(a), m_fes(fes) {
        init();
    }

    ~QoIFunctional () {}

    void init() const;
    double operator()(const GridFunction& f) const;
    double operator()(const Vector& f) const;

private:
    double m_pow = 1;
    FiniteElementSpace *m_fes = nullptr;
    mutable Vector m_bObs;
};

// Quantity of interest xt-functional; (x^a, fluxQ)
class QoIXtFunctional
{
public:

    QoIXtFunctional(FiniteElementSpace *tFes, FiniteElementSpace *xFes)
        : m_tFes(tFes), m_xFes(xFes) {
        init();
    }

    QoIXtFunctional(FiniteElementSpace *tFes, FiniteElementSpace *xFes,
                    double a)
        : m_pow(a), m_tFes(tFes), m_xFes(xFes) {
        init();
    }

    ~QoIXtFunctional () {}

    void init() const;
    double operator()(const Vector& f) const;

private:
    void add_vector(const Array<int> &vdofs,
                    const Vector& elvec,
                    Vector& b) const;

    double m_pow = 1;
    FiniteElementSpace *m_tFes = nullptr;
    FiniteElementSpace *m_xFes = nullptr;
    mutable Vector m_bObs;
};

}

#endif // HEAT_FUNCTIONALS_HPP
