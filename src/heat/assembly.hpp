#ifndef HEAT_ASSEMBLY_HPP
#define HEAT_ASSEMBLY_HPP

#include "mfem.hpp"


namespace heat {

/**
 * @brief Gradient Integrator in time; (du/dt, v)
 */
class TemporalGradientIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    TemporalGradientIntegrator() {}

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
class SpatialStiffnessIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    SpatialStiffnessIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    SpatialStiffnessIntegrator(mfem::MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    SpatialStiffnessIntegrator(mfem::MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with scalar material coefficient passed as a pointer
    SpatialStiffnessIntegrator(mfem::Coefficient *q_)
        : q(q_) {}

    //! Constructor with scalar material coefficient passed by reference
    SpatialStiffnessIntegrator(mfem::Coefficient& q_)
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
 * @brief VectorFE Stiffness Integrator in space; (div(u), div(v))
 * Uses, for example, Raviart-Thomas spaces.
 */
class SpatialVectorFEStiffnessIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    SpatialVectorFEStiffnessIntegrator() {}

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
 * @brief VectorFE Gradient Integrator in space; (grad(u), v)
 * Uses, for example, Raviart-Thomas spaces.
 */
class SpatialVectorFEGradientIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    SpatialVectorFEGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    SpatialVectorFEGradientIntegrator(mfem::MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    SpatialVectorFEGradientIntegrator(mfem::MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with material coefficent scalar passed as a pointer
    SpatialVectorFEGradientIntegrator(mfem::Coefficient *q_)
        : q(q_) {}

    //! Constructor with material coefficent scalar passed by reference
    SpatialVectorFEGradientIntegrator(mfem::Coefficient& q_)
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
 * @brief VectorFE Divergence LF Integrator in space; (f, div(v))
 * Uses, for example, Raviart-Thomas spaces.
 */
class SpatialVectorFEDivergenceLFIntegrator
        : public mfem::LinearFormIntegrator
{
public:
    //! Constructor with a coefficient and a scalar passed as arguments
    SpatialVectorFEDivergenceLFIntegrator(mfem::Coefficient& f, double a)
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
 * @brief Vector Stiffness Integrator in space; (div(u), div(v))
 */
class SpatialVectorStiffnessIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    SpatialVectorStiffnessIntegrator() {}

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
 * @brief Vector Gradient Integrator in space; (grad(u), v)
 */
class SpatialVectorGradientIntegrator
        : public mfem::BilinearFormIntegrator
{
public:
    SpatialVectorGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    SpatialVectorGradientIntegrator(mfem::MatrixCoefficient *M_)
        : M(M_) {}

    //! Constructor with material coefficent matrix passed by reference
    SpatialVectorGradientIntegrator(mfem::MatrixCoefficient& M_)
        : M(&M_) {}

    //! Constructor with material coefficent scalar passed as a pointer
    SpatialVectorGradientIntegrator(mfem::Coefficient *q_)
        : q(q_) {}

    //! Constructor with material coefficent scalar passed by reference
    SpatialVectorGradientIntegrator(mfem::Coefficient& q_)
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
 * @brief Vector Divergence LF Integrator in space; (f, div(v))
 */
class SpatialVectorDivergenceLFIntegrator
        : public mfem::LinearFormIntegrator
{
public:
    //! Constructor with a coefficient and a scalar passed as arguments
    SpatialVectorDivergenceLFIntegrator(mfem::Coefficient& f, double a)
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
