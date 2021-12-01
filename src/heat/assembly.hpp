#ifndef HEAT_ASSEMBLY_HPP
#define HEAT_ASSEMBLY_HPP

#include "mfem.hpp"


namespace heat {

/**
 * @brief Gradient Integrator in time; (du/dt, v)
 */
class GradientIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    GradientIntegrator() {}

    //! Assembles the Gradient integrator on a given temporal mesh element
    void AssembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    void AssembleElementMatrix2
    (const mfem::FiniteElement &,
     const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
#endif
};

/**
 * @brief Stiffness Integrator in space; (Q grad(u), Q grad(v))
 */
class StiffnessIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    StiffnessIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    StiffnessIntegrator(mfem::MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    StiffnessIntegrator(mfem::MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with scalar material coefficient passed as a pointer
    StiffnessIntegrator(mfem::Coefficient *q_)
        : q(q_) {}

    //! Constructor with scalar material coefficient passed by reference
    StiffnessIntegrator(mfem::Coefficient& q_)
        : q(&q_) {}

    //! Assembles the Stiffness integrator on a given spatial mesh element
    void AssembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    void AssembleElementMatrix2
    (const mfem::FiniteElement &,
     const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &)
    override {}

private:
    mfem::MatrixCoefficient *M = nullptr;
    mfem::Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    mfem::DenseMatrix dshape;
#endif
};

/**
 * @brief VectorFE Stiffness Integrator; (div(u), div(v))
 * Uses, for example, Raviart-Thomas spaces.
 */
class VectorFEStiffnessIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    VectorFEStiffnessIntegrator() {}

    //! Assembles the VectorFE Stiffness integrator on a given spatial mesh element
    void AssembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    void AssembleElementMatrix2
    (const mfem::FiniteElement &,
     const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    mfem::Vector divshape;
#endif
};

/**
 * @brief VectorFE Gradient Integrator; (grad(u), v)
 * Uses, for example, Raviart-Thomas spaces.
 */
class VectorFEGradientIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    VectorFEGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    VectorFEGradientIntegrator(mfem::MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    VectorFEGradientIntegrator(mfem::MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with material coefficent scalar passed as a pointer
    VectorFEGradientIntegrator(mfem::Coefficient *q_)
        : q(q_) {}

    //! Constructor with material coefficent scalar passed by reference
    VectorFEGradientIntegrator(mfem::Coefficient& q_)
        : q(&q_) {}

    void AssembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override {}

    //! Assembles the VectorFE Gradient integrator on a given spatial mesh element
    //! The test and trial funtions are different
    void AssembleElementMatrix2
    (const mfem::FiniteElement &,
     const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &)
    override;

private:
    mfem::MatrixCoefficient *M = nullptr;
    mfem::Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    mfem::DenseMatrix gradshape, vshape;
#endif
};

/**
 * @brief VectorFE Divergence LF Integrator; (f, div(v))
 * Uses, for example, Raviart-Thomas spaces.
 */
class VectorFEDivergenceLFIntegrator
        : public mfem::LinearFormIntegrator
{
public:
    //! Constructor with a coefficient and a scalar passed as arguments
    VectorFEDivergenceLFIntegrator(mfem::Coefficient& f, double a)
        : F(f), m_coeff(a) {}

    //! Assembles the VectorFE linear form integrator on a given spatial mesh element
    void AssembleRHSElementVect (const mfem::FiniteElement &,
                                 mfem::ElementTransformation &,
                                 mfem::Vector &) override;

private:
    mfem::Coefficient &F;
    double m_coeff;
#ifndef MFEM_THREAD_SAFE
    mfem::Vector divshape;
#endif
};

/**
 * @brief mfem::Vector Stiffness Integrator; (div(u), div(v))
 */
class VectorStiffnessIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    VectorStiffnessIntegrator() {}

    //! Assembles the mfem::Vector Stiffness integrator on a given spatial mesh element
    void AssembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    void AssembleElementMatrix2
    (const mfem::FiniteElement &,
     const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    mfem::Vector divshape;
    mfem::DenseMatrix dshape;
#endif
};

/**
 * @brief mfem::Vector Gradient Integrator; (grad(u), v)
 */
class VectorGradientIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    VectorGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    VectorGradientIntegrator(mfem::MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    VectorGradientIntegrator(mfem::MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with material coefficent scalar passed as a pointer
    VectorGradientIntegrator(mfem::Coefficient *q_)
        : q(q_) {}

    //! Constructor with material coefficent scalar passed by reference
    VectorGradientIntegrator(mfem::Coefficient& q_)
        : q(&q_) {}

    void AssembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override {}

    //! Assembles the VectorFE Gradient integrator on a given spatial mesh element
    //! The test and trial funtions are different
    void AssembleElementMatrix2
    (const mfem::FiniteElement &,
     const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &)
    override;

private:
    void AssembleBlock(const int dim,
                       const int test_ndofs,
                       const int trial_ndofs,
                       const mfem::Vector& shape,
                       const mfem::DenseMatrix& dshape,
                       mfem::DenseMatrix& outMat) const;

    mfem::MatrixCoefficient *M = nullptr;
    mfem::Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
#endif
};

/**
 * @brief mfem::Vector Divergence LF Integrator; (f, div(v))
 */
class VectorDivergenceLFIntegrator
        : public mfem::LinearFormIntegrator
{
public:
    //! Constructor with a coefficient and a scalar passed as arguments
    VectorDivergenceLFIntegrator(mfem::Coefficient& f, double a)
        : F(f), m_coeff(a) {}

    //! Assembles the mfem::Vector linear form integrator on a given spatial mesh element
    void AssembleRHSElementVect (const mfem::FiniteElement &,
                                 mfem::ElementTransformation &,
                                 mfem::Vector &) override;

private:
    mfem::Coefficient &F;
    double m_coeff;
#ifndef MFEM_THREAD_SAFE
    mfem::Vector divshape;
    mfem::DenseMatrix dshape;
#endif
};

}

#endif // HEAT_ASSEMBLY_HPP
