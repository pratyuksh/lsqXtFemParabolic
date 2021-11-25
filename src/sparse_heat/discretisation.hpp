#ifndef SPARSE_HEAT_DISCRETISATION_HPP
#define SPARSE_HEAT_DISCRETISATION_HPP

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../core/config.hpp"
#include "../heat/test_cases.hpp"
#include "../mymfem/nested_hierarchy.hpp"


namespace sparseHeat {

/**
 * @brief Base class for Least-squares sparse space-time FE discretisation of the heat equation
 */
class LsqSparseXtFem
{
public:
    //! Constructor with config and test case as arguments
    LsqSparseXtFem (const nlohmann::json&,
                    std::shared_ptr<heat::TestCases>&);

    //! Constructor with config, test case
    //! and mesh hierarchy as arguments
    LsqSparseXtFem (const nlohmann::json&,
                    std::shared_ptr<heat::TestCases>&,
                    std::shared_ptr<mymfem::NestedMeshHierarchy>&);

    //! Default destructor
    virtual ~LsqSparseXtFem ();

    //! Sets nested hierarchies of FE spaces and spatial boundary dofs
    void setNestedFEHierarchyAndSpatialBoundaryDofs
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&);

    //! Sets nested hierarchies of FE spaces
    void setNestedFEHierarchy
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&);

    void setNestedFEHierarchyForTemperature
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&);

    virtual void setNestedFEHierarchyForHeatFlux
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&) = 0;

    //! Sets the spatial boundary dofs
    void setSpatialBoundaryDofs();

    //! Assembles the sub-matrices in the discretisation
    void assembleSystemSubMatrices();

protected:
    void assembleTemporalInitial();
    void assembleTemporalMass();
    void assembleTemporalStiffness();
    void assembleTemporalGradient();

    void assembleSpatialMassForTemperature();
    void assembleSpatialStiffnessForTemperature();

    virtual void assembleSpatialMassForHeatFlux() = 0;
    virtual void assembleSpatialStiffnessForHeatFlux() = 0;
    virtual void assembleSpatialGradient() = 0;
    virtual void assembleSpatialDivergence() = 0;

public:
    //! Builds the system matrix
    void buildSystemMatrix();

private:
    //! Builds the system matrix blocks
    void buildSystemBlocks();

    void buildSystemBlock11();
    void buildSystemBlock22();
    void buildSystemBlock12();
    void buildSystemBlock21();

    inline int getSpatialIndex(int i) const;

public:
    void assembleRhs(std::shared_ptr<BlockVector>&) const;

private:
    void assembleICs(Vector& b) const;

    void assembleSource(std::shared_ptr<BlockVector>& B) const;
    void assembleSourceWithTemporalGradientOfTemperatureBasis(Vector& b) const;
    void assembleSourceWithTemporalGradientOfTemperatureBasisAtGivenTime
    (double t,
     const std::shared_ptr<FiniteElementSpace>&,
     Vector&) const;

    void assembleSourceWithSpatialDivergenceOfHeatFluxBasis(Vector& b) const;

protected:
    virtual void assembleSourceWithSpatialDivergenceOfHeatFluxBasisAtGivenTime
    (double t,
     const std::shared_ptr<FiniteElementSpace>&,
     Vector&) const = 0;

public:
    void applyBCs(SparseMatrix& A) const;

    void applyBCs(BlockVector& B) const;

public:
    SparseMatrix* getSystemMatrix() const {
        return m_systemMatrix;
    }

    SparseMatrix* getSystemBlock11() const {
        return m_systemBlock11;
    }

    SparseMatrix* getSystemBlock12() const {
        return m_systemBlock12;
    }

    SparseMatrix* getSystemBlock21() const {
        return m_systemBlock21;
    }

    SparseMatrix* getSystemBlock22() const {
        return m_systemBlock22;
    }

    std::shared_ptr<mymfem::NestedFEHierarchy>
    getSpatialNestedFEHierarchyForTemperature() const {
        return m_spatialNestedFEHierarchyTemperature;
    }

    std::shared_ptr<mymfem::NestedFEHierarchy>
    getSpatialNestedFEHierarchyForHeatFlux() const {
        return m_spatialNestedFEHierarchyHeatFlux;
    }

protected:
    const nlohmann::json& m_config;
    std::shared_ptr<heat::TestCases> m_testCase;

    int m_deg;

    double m_endTime;
    int m_minTemporalLevel, m_maxTemporalLevel, m_numLevels;

    int m_xDim;
    Array<FiniteElementCollection*> m_xFecs;
    std::shared_ptr<mymfem::NestedFEHierarchy> m_spatialNestedFEHierarchyTemperature;
    std::shared_ptr<mymfem::NestedFEHierarchy> m_spatialNestedFEHierarchyHeatFlux;

    Array<int> m_essDofs;

    std::shared_ptr<BlockMatrix> m_temporalMass;
    std::shared_ptr<BlockMatrix> m_temporalStiffness;
    std::shared_ptr<BlockMatrix> m_temporalGradient;
    std::shared_ptr<BlockMatrix> m_temporalInitial;

    std::shared_ptr<BlockMatrix> m_spatialMass1, m_spatialMass2;
    std::shared_ptr<BlockMatrix> m_spatialStiffness1, m_spatialStiffness2;
    std::shared_ptr<BlockMatrix> m_spatialGradient;
    std::shared_ptr<BlockMatrix> m_spatialDivergence;

    SparseMatrix* m_systemBlock11 = nullptr;
    SparseMatrix* m_systemBlock12 = nullptr;
    SparseMatrix* m_systemBlock21 = nullptr;
    SparseMatrix* m_systemBlock22 = nullptr;

    SparseMatrix* m_systemMatrix = nullptr;
};


/**
 * @brief Least-squares sparse space-time FEM for the heat equation,
 * with H1 discretisation in time and
 * Hdiv discretisation in space
 */
class LsqSparseXtFemH1Hdiv : public LsqSparseXtFem
{
public:
    LsqSparseXtFemH1Hdiv (const nlohmann::json& config,
                          std::shared_ptr<heat::TestCases>& testCase);

    LsqSparseXtFemH1Hdiv
    (const nlohmann::json& config,
     std::shared_ptr<heat::TestCases>& testCase,
     std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialMeshHierarchy);

    void setNestedFEHierarchyForHeatFlux
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&) override;

    void assembleSpatialMassForHeatFlux() override;
    void assembleSpatialStiffnessForHeatFlux() override;
    void assembleSpatialGradient() override;
    void assembleSpatialDivergence() override;

private:
    void assembleSourceWithSpatialDivergenceOfHeatFluxBasisAtGivenTime
    (double t,
     const std::shared_ptr<FiniteElementSpace>&,
     Vector&) const override;
};

}

#endif // SPARSE_HEAT_DISCRETISATION_HPP
