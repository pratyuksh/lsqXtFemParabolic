#ifndef HEAT_ASSEMBLY_HPP
#define HEAT_ASSEMBLY_HPP

#include "mfem.hpp"
using namespace mfem;

namespace heat {

/**
 * @brief Gradient Integrator in time; (du/dt, v)
 */
class GradientIntegrator
        : public BilinearFormIntegrator
{
public:
    GradientIntegrator() {}

    //! Assembles the Gradient integrator on a given temporal mesh element
    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    Vector shape;
    DenseMatrix dshape;
#endif
};

/**
 * @brief Stiffness Integrator in space; (Q grad(u), Q grad(v))
 */
class StiffnessIntegrator
        : public BilinearFormIntegrator
{
public:
    StiffnessIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    StiffnessIntegrator(MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    StiffnessIntegrator(MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with scalar material coefficient passed as a pointer
    StiffnessIntegrator(Coefficient *q_)
        : q(q_) {}

    //! Constructor with scalar material coefficient passed by reference
    StiffnessIntegrator(Coefficient& q_)
        : q(&q_) {}

    //! Assembles the Stiffness integrator on a given spatial mesh element
    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
    MatrixCoefficient *M = nullptr;
    Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    DenseMatrix dshape;
#endif
};

/**
 * @brief VectorFE Stiffness Integrator; (div(u), div(v))
 * Uses, for example, Raviart-Thomas spaces.
 */
class VectorFEStiffnessIntegrator
        : public BilinearFormIntegrator
{
public:
    VectorFEStiffnessIntegrator() {}

    //! Assembles the VectorFE Stiffness integrator on a given spatial mesh element
    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    Vector divshape;
#endif
};

/**
 * @brief VectorFE Gradient Integrator; (grad(u), v)
 * Uses, for example, Raviart-Thomas spaces.
 */
class VectorFEGradientIntegrator
        : public BilinearFormIntegrator
{
public:
    VectorFEGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    VectorFEGradientIntegrator(MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    VectorFEGradientIntegrator(MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with material coefficent scalar passed as a pointer
    VectorFEGradientIntegrator(Coefficient *q_)
        : q(q_) {}

    //! Constructor with material coefficent scalar passed by reference
    VectorFEGradientIntegrator(Coefficient& q_)
        : q(&q_) {}

    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override {}

    //! Assembles the VectorFE Gradient integrator on a given spatial mesh element
    //! The test and trial funtions are different
    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override;

private:
    MatrixCoefficient *M = nullptr;
    Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    DenseMatrix gradshape, vshape;
#endif
};

/**
 * @brief VectorFE Divergence LF Integrator; (f, div(v))
 * Uses, for example, Raviart-Thomas spaces.
 */
class VectorFEDivergenceLFIntegrator
        : public LinearFormIntegrator
{
public:
    //! Constructor with a coefficient and a scalar passed as arguments
    VectorFEDivergenceLFIntegrator(Coefficient& f, double a)
        : F(f), m_coeff(a) {}

    //! Assembles the VectorFE linear form integrator on a given spatial mesh element
    void AssembleRHSElementVect (const FiniteElement &,
                                 ElementTransformation &,
                                 Vector &) override;

private:
    Coefficient &F;
    double m_coeff;
#ifndef MFEM_THREAD_SAFE
    Vector divshape;
#endif
};

/**
 * @brief Vector Stiffness Integrator; (div(u), div(v))
 */
class VectorStiffnessIntegrator
        : public BilinearFormIntegrator
{
public:
    VectorStiffnessIntegrator() {}

    //! Assembles the Vector Stiffness integrator on a given spatial mesh element
    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    Vector divshape;
    DenseMatrix dshape;
#endif
};

/**
 * @brief Vector Gradient Integrator; (grad(u), v)
 */
class VectorGradientIntegrator
        : public BilinearFormIntegrator
{
public:
    VectorGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    VectorGradientIntegrator(MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    VectorGradientIntegrator(MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with material coefficent scalar passed as a pointer
    VectorGradientIntegrator(Coefficient *q_)
        : q(q_) {}

    //! Constructor with material coefficent scalar passed by reference
    VectorGradientIntegrator(Coefficient& q_)
        : q(&q_) {}

    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override {}

    //! Assembles the VectorFE Gradient integrator on a given spatial mesh element
    //! The test and trial funtions are different
    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override;

private:
    void AssembleBlock(const int dim,
                       const int test_ndofs,
                       const int trial_ndofs,
                       const Vector& shape,
                       const DenseMatrix& dshape,
                       DenseMatrix& outMat) const;

    MatrixCoefficient *M = nullptr;
    Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    Vector shape;
    DenseMatrix dshape;
#endif
};

/**
 * @brief Vector Divergence LF Integrator; (f, div(v))
 */
class VectorDivergenceLFIntegrator
        : public LinearFormIntegrator
{
public:
    //! Constructor with a coefficient and a scalar passed as arguments
    VectorDivergenceLFIntegrator(Coefficient& f, double a)
        : F(f), m_coeff(a) {}

    //! Assembles the Vector linear form integrator on a given spatial mesh element
    void AssembleRHSElementVect (const FiniteElement &,
                                 ElementTransformation &,
                                 Vector &) override;

private:
    Coefficient &F;
    double m_coeff;
#ifndef MFEM_THREAD_SAFE
    Vector divshape;
    DenseMatrix dshape;
#endif
};

}

#endif // HEAT_ASSEMBLY_HPP
