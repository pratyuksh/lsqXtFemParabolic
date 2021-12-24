#ifndef HEAT_DISCRETISATION_HPP
#define HEAT_DISCRETISATION_HPP

#include "mfem.hpp"

#include <iostream>

#include "../core/config.hpp"
#include "test_cases.hpp"


namespace heat {

/**
 * @brief Base class for Least-squares space-time FE discretisation
 * of the heat equation
 */
class LsqXtFem
{
public:
    LsqXtFem (const nlohmann::json&, std::shared_ptr<heat::TestCases>&);

    ~LsqXtFem ();

    //! Releases all the memory allocated by this class
    void reset(bool firstPass=true);

    void resetFeSpaces();
    void resetSystemOperators();

    void resetSystemBlocks();
    void resetMediumIndependentSystemBlocks();
    void resetMediumDependentSystemBlocks();

    void resetSystemSubMatrices();
    void resetMediumIndependentSystemSubMatrices();
    void resetMediumDependentSystemSubMatrices();
    
    //! Sets the FE spaces, block offsets and boundary conditions
    void resetFeSpacesAndBlockOffsetsAndSpatialBoundaryDofs
        (std::shared_ptr<mfem::Mesh>&, std::shared_ptr<mfem::Mesh>&);

    void setFeSpacesAndBlockOffsetsAndSpatialBoundaryDofs
        (std::shared_ptr<mfem::Mesh>&, std::shared_ptr<mfem::Mesh>&);

    //! Sets temporal and spatial FE spaces for temperature
    //! and heat flux
    void setFeSpaces
    (std::shared_ptr<mfem::Mesh>&,
     std::shared_ptr<mfem::Mesh>&);

    void setSpatialFeSpaceForTemperature
    (std::shared_ptr<mfem::Mesh>&);

    virtual void setSpatialFeSpaceForHeatFlux
    (std::shared_ptr<mfem::Mesh>&) = 0;

    //! Sets the block offsets for the block-structured system
    void setBlockOffsets();

    void setSpatialBoundaryDofs();

    //! Assembles the sub-matrices in the discretisation
    void reassembleSystemSubMatrices();
    void assembleSystemSubMatrices();
    void assembleMaterialIndependentSystemSubMatrices();
    void assembleMaterialDependentSystemSubMatrices();

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

    //! Builds the linear system matrix from the blocks
    //! as an MFEM operator
    void rebuildSystemOperator();
    void buildSystemOperator();

    //! Builds the linear system matrix from the blocks
    //! as a monolithic matrix
    void rebuildSystemMatrix();
    void buildSystemMatrix();

    //! Builds the upper-triangle linear system matrix
    //! from the blocks as a monolithic matrix
    void rebuildUpperTriangleOfSystemMatrix();
    void buildUpperTriangleOfSystemMatrix();

private:
    //! Builds the system matrix blocks
    void rebuildSystemBlocks();
    void buildSystemBlocks();

    void buildSystemBlock22();

    void rebuildSystemBlock11();
    void buildSystemBlock11();
    void buildMediumIndependentSystemBlock11();
    void buildMediumDependentSystemBlock11();

    void rebuildSystemBlock12AndBlock21();
    void buildSystemBlock12AndBlock21();
    void buildMediumIndependentSystemBlock12();
    void buildMediumDependentSystemBlock12();

public:
    void assembleRhs(mfem::BlockVector *) const;

protected:
    void assembleICs(mfem::Vector&) const;
    void assembleSpatialICs(mfem::Vector&) const;

    void assembleSource(mfem::BlockVector&) const;
    void assembleSourceWithTemporalGradientOfTemperatureBasis
    (mfem::Vector&) const;
    void assembleSpatialSourceWithTemperatureBasis
    (mfem::Vector&, double) const;
    void assembleSourceWithSpatialDivergenceOfHeatFluxBasis(
            mfem::Vector&) const;

    virtual void assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis
    (mfem::Vector&, double) const = 0;

    void addVector(const mfem::Array<int> &,
                   const mfem::Vector&, mfem::Vector&) const;

    //! Applies boundary conditions to a SparseMatrix
    void applyBCs(mfem::SparseMatrix&) const;

    //! Applies boundary conditions to a BlockVector
    void applyBCs(mfem::BlockVector&) const;

public:
    std::shared_ptr<heat::TestCases> getTestCase() const {
        return m_testCase;
    }
    
    mfem::Mesh* getTemporalMesh() const {
        return m_temporalFeSpace->GetMesh();
    }
    
    mfem::Mesh* getSpatialMesh() const {
        return m_spatialFeSpaces[0]->GetMesh();
    }

    mfem::FiniteElementCollection* getTemporalFeCollection() const {
        return m_temporalFeCollection;
    }

    mfem::Array<mfem::FiniteElementCollection*>
    getSpatialFeCollections() const {
        return m_spatialFeCollections;
    }
    
    mfem::FiniteElementSpace* getTemporalFeSpace() const {
        return m_temporalFeSpace;
    }
    
    mfem::Array<mfem::FiniteElementSpace*> getSpatialFeSpaces() const {
        return m_spatialFeSpaces;
    }

    mfem::BlockOperator* getSystemOperator() const {
        return m_systemOperator;
    }

    mfem::Array<int> getBlockOffsets() const {
        return m_blockOffsets;
    }

    mfem::SparseMatrix* getSystemMatrix() const {
        return m_systemMatrix;
    }
    
    //! Prints some info
    void print() const {
        std::cout << "Degree in x and t: "
                  << m_deg << std::endl;
    }

    mfem::SparseMatrix* getSystemBlock11() const {
        return m_systemBlock11;
    }

    mfem::SparseMatrix* getSystemBlock12() const {
        return m_systemBlock12;
    }

    mfem::SparseMatrix* getSystemBlock21() const {
        return m_systemBlock21;
    }

    mfem::SparseMatrix* getSystemBlock22() const {
        return m_systemBlock22;
    }

protected:
    const nlohmann::json& m_config;
    
    int m_xDim, m_deg;
    std::shared_ptr<heat::TestCases> m_testCase;
    
    mfem::FiniteElementCollection* m_temporalFeCollection = nullptr;
    mfem::FiniteElementSpace* m_temporalFeSpace = nullptr;
    
    mfem::Array<mfem::FiniteElementCollection*> m_spatialFeCollections;
    mfem::Array<mfem::FiniteElementSpace*> m_spatialFeSpaces;
    mfem::Array<int> m_blockOffsets;
    
    mfem::SparseMatrix *m_temporalMass = nullptr;
    mfem::SparseMatrix *m_temporalStiffness = nullptr;
    mfem::SparseMatrix *m_temporalGradient = nullptr;
    mfem::SparseMatrix *m_temporalInitial = nullptr;

    mfem::SparseMatrix *m_spatialMass1 = nullptr;
    mfem::SparseMatrix *m_spatialMass2 = nullptr;
    mfem::SparseMatrix *m_spatialStiffness1 = nullptr;
    mfem::SparseMatrix *m_spatialStiffness2 = nullptr;
    mfem::SparseMatrix *m_spatialGradient = nullptr;
    mfem::SparseMatrix *m_spatialDivergence = nullptr;

    mfem::SparseMatrix *m_systemBlock11 = nullptr;
    mfem::SparseMatrix *m_systemBlock12 = nullptr;
    mfem::SparseMatrix *m_systemBlock21 = nullptr;
    mfem::SparseMatrix *m_systemBlock22 = nullptr;

    mfem::BlockOperator *m_systemOperator = nullptr;
    mfem::SparseMatrix *m_systemMatrix = nullptr;

    mfem::Array<int> m_spatialEssentialBoundaryMarker;
    mfem::Array<int> m_spatialEssentialDofs;
    mfem::Array<int> m_essentialDofs;

    bool m_firstPassOfAssembleSystemSubMatrices = true;
    bool m_firstPassOfAssembleSystemBlocks = true;
    mfem::SparseMatrix *m_materialIndependentSystemBlock11 = nullptr;
    mfem::SparseMatrix *m_materialIndependentSystemBlock12 = nullptr;
};

/**
 * @brief Least-squares space-time FEM for the heat equation,
 * with H1 discretisation in both time and space
 */
class LsqXtFemH1H1 : public LsqXtFem
{
public:
    LsqXtFemH1H1 (const nlohmann::json& config,
                  std::shared_ptr<heat::TestCases>& testCase)
        : LsqXtFem(config, testCase) {}

    void setSpatialFeSpaceForHeatFlux
        (std::shared_ptr<mfem::Mesh>&) override;

    void assembleSpatialMassForHeatFlux() override;
    void assembleSpatialStiffnessForHeatFlux() override;
    void assembleSpatialGradient() override;
    void assembleSpatialDivergence() override;

private:
    void assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis
    (mfem::Vector&, double) const override;
};

/**
 * @brief Least-squares space-time FEM for the heat equation,
 * with H1 discretisation in time and
 * Hdiv discretisation in space
 */class LsqXtFemH1Hdiv : public LsqXtFem
{
public:
    LsqXtFemH1Hdiv (const nlohmann::json& config,
                        std::shared_ptr<heat::TestCases>& testCase)
        : LsqXtFem(config, testCase) {}

    void setSpatialFeSpaceForHeatFlux
        (std::shared_ptr<mfem::Mesh>&) override;

    void assembleSpatialMassForHeatFlux() override;
    void assembleSpatialStiffnessForHeatFlux() override;
    void assembleSpatialGradient() override;
    void assembleSpatialDivergence() override;

private:
    void assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis
    (mfem::Vector&, double) const override;
};

}

#endif // HEAT_DISCRETISATION_HPP
