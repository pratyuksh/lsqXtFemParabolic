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

    //! Default destructor
    virtual ~LsqSparseXtFem ();

    //! Sets nested hierarchies of FE spaces
    void setNestedFEHierarchy
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&);

    void setNestedFEHierarchyForTemperature
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&);

    virtual void setNestedFEHierarchyForHeatFlux
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&) = 0;

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
    //! Builds the system matrix blocks
    void buildSystemBlocks();

private:
    void buildSystemBlock11();
    void buildSystemBlock22();
    void buildSystemBlock12();
    void buildSystemBlock21();

    inline int getSpatialIndex(int i) {
        return m_numLevels - i - 1;
    }

public:
    SparseMatrix* getSystemBlock11() {
        return m_block11;
    }

    SparseMatrix* getSystemBlock12() {
        return m_block12;
    }

    SparseMatrix* getSystemBlock21() {
        return m_block21;
    }

    SparseMatrix* getSystemBlock22() {
        return m_block22;
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

    std::shared_ptr<BlockMatrix> m_temporalMass;
    std::shared_ptr<BlockMatrix> m_temporalStiffness;
    std::shared_ptr<BlockMatrix> m_temporalGradient;
    std::shared_ptr<BlockMatrix> m_temporalInitial;

    std::shared_ptr<BlockMatrix> m_spatialMass1, m_spatialMass2;
    std::shared_ptr<BlockMatrix> m_spatialStiffness1, m_spatialStiffness2;
    std::shared_ptr<BlockMatrix> m_spatialGradient;
    std::shared_ptr<BlockMatrix> m_spatialDivergence;

    SparseMatrix* m_block11 = nullptr;
    SparseMatrix* m_block12 = nullptr;
    SparseMatrix* m_block21 = nullptr;
    SparseMatrix* m_block22 = nullptr;
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
                          std::shared_ptr<heat::TestCases>& testCase)
        : LsqSparseXtFem(config, testCase) {}

    void setNestedFEHierarchyForHeatFlux
    (std::shared_ptr<mymfem::NestedMeshHierarchy>&) override;

    void assembleSpatialMassForHeatFlux() override;
    void assembleSpatialStiffnessForHeatFlux() override;
    void assembleSpatialGradient() override;
    void assembleSpatialDivergence() override;
};

}

#endif // SPARSE_HEAT_DISCRETISATION_HPP
