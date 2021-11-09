#ifndef SPARSE_HEAT_ASSEMBLY_HPP
#define SPARSE_HEAT_ASSEMBLY_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../mymfem/my_bilinearForms.hpp"


namespace sparseHeat {

/**
 * @brief Mass Integrator; (u, v)
 */
class MassIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    //! Default constructor
    explicit MassIntegrator() {}

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
class StiffnessIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    explicit StiffnessIntegrator() {}

    //! Constructor with material coefficent matrix passed as a pointer
    explicit StiffnessIntegrator(MatrixCoefficient *M)
        : m_matrixCoeff(M) {}

    //! Constructor with material coefficent matrix passed by reference
    explicit StiffnessIntegrator(MatrixCoefficient& M)
        : m_matrixCoeff(&M) {}

    //! Constructor with scalar material coefficient passed as a pointer
    explicit StiffnessIntegrator(Coefficient *q)
        : m_scalarCoeff(q) {}

    //! Constructor with scalar material coefficient passed by reference
    explicit StiffnessIntegrator(Coefficient& q)
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
class VectorFEMassIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    //! Default constructor
    explicit VectorFEMassIntegrator() {}

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
class VectorFEStiffnessIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    //! Default constructor
    explicit VectorFEStiffnessIntegrator() {}

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

}

#endif // SPARSE_HEAT_ASSEMBLY_HPP
