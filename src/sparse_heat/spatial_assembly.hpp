#ifndef SPARSE_HEAT_SPATIAL_ASSEMBLY_HPP
#define SPARSE_HEAT_SPATIAL_ASSEMBLY_HPP

#include "mfem.hpp"
using namespace mfem;

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
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    //! Assembles the mass integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;
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
    explicit SpatialStiffnessIntegrator(MatrixCoefficient *M)
        : m_matrixCoeff(M) {}

    //! Constructor with material coefficent matrix passed by reference
    explicit SpatialStiffnessIntegrator(MatrixCoefficient& M)
        : m_matrixCoeff(&M) {}

    //! Constructor with scalar material coefficient passed as a pointer
    explicit SpatialStiffnessIntegrator(Coefficient *q)
        : m_scalarCoeff(q) {}

    //! Constructor with scalar material coefficient passed by reference
    explicit SpatialStiffnessIntegrator(Coefficient& q)
        : m_scalarCoeff(&q) {}

    //! Assembles the stiffness integrator on an element
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    //! Assembles the stiffness integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

private:
    MatrixCoefficient *m_matrixCoeff = nullptr;
    Coefficient *m_scalarCoeff = nullptr;
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
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    //! Assembles the vectorFE mass integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;
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
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    //! Assembles the vectorFE stiffness integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;
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
    explicit SpatialVectorFEGradientIntegrator(MatrixCoefficient *M)
        : m_matrixCoeff(M) {}

    //! Constructor with material coefficent matrix passed by reference
    explicit SpatialVectorFEGradientIntegrator(MatrixCoefficient& M)
        : m_matrixCoeff(&M) {}

    //! Constructor with scalar material coefficient passed as a pointer
    explicit SpatialVectorFEGradientIntegrator(Coefficient *q)
        : m_scalarCoeff(q) {}

    //! Constructor with scalar material coefficient passed by reference
    explicit SpatialVectorFEGradientIntegrator(Coefficient& q)
        : m_scalarCoeff(&q) {}

    //! Assembles the mixed scalar-vectorFE gradient integrator on an element
    void assembleElementMatrix
    (const FiniteElement &, const FiniteElement &,
     ElementTransformation &, DenseMatrix &) override;

    //! Assembles the mixed scalar-vectorFE gradient integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void assembleElementMatrix2
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

private:
    MatrixCoefficient *m_matrixCoeff = nullptr;
    Coefficient *m_scalarCoeff = nullptr;
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
    (const FiniteElement &, const FiniteElement &,
     ElementTransformation &, DenseMatrix &) override;

    //! Assembles the mixed scalar-vectorFE divergence integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void assembleElementMatrix2
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;
};

}

#endif // SPARSE_HEAT_ASSEMBLY_HPP
