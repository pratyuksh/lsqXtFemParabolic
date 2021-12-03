#ifndef SPARSE_HEAT_SPATIAL_ASSEMBLY_HPP
#define SPARSE_HEAT_SPATIAL_ASSEMBLY_HPP

#include "mfem.hpp"

#include "../mymfem/my_bilinearForms.hpp"


namespace sparseHeat {

/**
 * @brief Mass Integrator; (u, v)
 */
class SpatialMassIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit SpatialMassIntegrator() {}

    //! Assembles the mass integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    //! Assembles the mass integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;
};


/**
 * @brief Stiffness Integrator; (Q grad(u), Q grad(v))
 */
class SpatialStiffnessIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit SpatialStiffnessIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    explicit SpatialStiffnessIntegrator(mfem::MatrixCoefficient *M)
        : m_matrixCoeff(M) {}

    //! Constructor with material coefficent matrix passed by reference
    explicit SpatialStiffnessIntegrator(mfem::MatrixCoefficient& M)
        : m_matrixCoeff(&M) {}

    //! Constructor with scalar material coefficient passed as a pointer
    explicit SpatialStiffnessIntegrator(mfem::Coefficient *q)
        : m_scalarCoeff(q) {}

    //! Constructor with scalar material coefficient passed by reference
    explicit SpatialStiffnessIntegrator(mfem::Coefficient& q)
        : m_scalarCoeff(&q) {}

    //! Assembles the stiffness integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    //! Assembles the stiffness integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

private:
    mfem::MatrixCoefficient *m_matrixCoeff = nullptr;
    mfem::Coefficient *m_scalarCoeff = nullptr;
};


/**
 * @brief VectorFE Mass Integrator; (u, v)
 */
class SpatialVectorFEMassIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit SpatialVectorFEMassIntegrator() {}

    //! Assembles the vectorFE mass integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    //! Assembles the vectorFE mass integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;
};


/**
 * @brief VectorFE Stiffness Integrator; (div(u), div(v))
 */
class SpatialVectorFEStiffnessIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit SpatialVectorFEStiffnessIntegrator() {}

    //! Assembles the vectorFE stiffness integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    //! Assembles the vectorFE stiffness integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;
};


/**
 * @brief VectorFE Gradient Integrator; (grad(u), v)
 */
class SpatialVectorFEGradientIntegrator
        : public mymfem::BlockMixedBilinearFormIntegrator
{
public:
    explicit SpatialVectorFEGradientIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    explicit SpatialVectorFEGradientIntegrator(mfem::MatrixCoefficient *M)
        : m_matrixCoeff(M) {}

    //! Constructor with material coefficent matrix passed by reference
    explicit SpatialVectorFEGradientIntegrator(mfem::MatrixCoefficient& M)
        : m_matrixCoeff(&M) {}

    //! Constructor with scalar material coefficient passed as a pointer
    explicit SpatialVectorFEGradientIntegrator(mfem::Coefficient *q)
        : m_scalarCoeff(q) {}

    //! Constructor with scalar material coefficient passed by reference
    explicit SpatialVectorFEGradientIntegrator(mfem::Coefficient& q)
        : m_scalarCoeff(&q) {}

    //! Assembles the mixed scalar-vectorFE gradient integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &) override;

    //! Assembles the mixed scalar-vectorFE gradient integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    void assembleElementMatrix2
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

private:
    mfem::MatrixCoefficient *m_matrixCoeff = nullptr;
    mfem::Coefficient *m_scalarCoeff = nullptr;
};


/**
 * @brief VectorFE Divergence Integrator; (div(u), v)
 */
class SpatialVectorFEDivergenceIntegrator
        : public mymfem::BlockMixedBilinearFormIntegrator
{
public:
    explicit SpatialVectorFEDivergenceIntegrator() {}

    //! Assembles the mixed scalar-vectorFE divergence integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, const mfem::FiniteElement &,
     mfem::ElementTransformation &, mfem::DenseMatrix &) override;

    //! Assembles the mixed scalar-vectorFE divergence integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    void assembleElementMatrix2
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;
};


/**
 * @brief Vector Mass Integrator; (u, v)
 */
class SpatialVectorMassIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit SpatialVectorMassIntegrator() {}

    //! Assembles the vector mass integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    //! Assembles the vector mass integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;
};


/**
 * @brief Vector Stiffness Integrator; (div(u), div(v))
 */
class SpatialVectorStiffnessIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit SpatialVectorStiffnessIntegrator() {}

    //! Assembles the vectorFE stiffness integrator on an element
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;

    //! Assembles the vectorFE stiffness integrator between two nested mesh elements
    void assembleElementMatrix
    (const mfem::FiniteElement &, mfem::ElementTransformation &,
     const mfem::FiniteElement &, mfem::ElementTransformation &,
     mfem::DenseMatrix &) override;
};

}

#endif // SPARSE_HEAT_ASSEMBLY_HPP
