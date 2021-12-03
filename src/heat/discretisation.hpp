#ifndef HEAT_DISCRETISATION_HPP
#define HEAT_DISCRETISATION_HPP

#include "mfem.hpp"

#include <iostream>

#include "../core/config.hpp"
#include "test_cases.hpp"


namespace heat {

/**
 * @brief Base class for Least-squares space-time FE discretisation of the heat equation
 */
class LsqXtFEM
{
public:
    //! Constructor with config and test case as arguments
    LsqXtFEM (const nlohmann::json&, std::shared_ptr<heat::TestCases>&);

    //! Default destructor
    ~LsqXtFEM ();
    
    //! Sets the FE space and boundary conditions
    virtual void set(std::shared_ptr<mfem::Mesh>&,
                     std::shared_ptr<mfem::Mesh>&) {
        std::cout << "Function 'set' not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    //! Assembles the operators in the discretisation
    //! that are independent of the material properties
    virtual void assembleSystemMediumIndependent()  {
        std::cout << "Function "
                     "'assemble_system_medium_independent' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    //! Assembles the operators in the discretisation
    //! that depend on the material properties
    virtual void assembleSystemMediumDependent()  {
        std::cout << "Function "
                     "'assemble_system_medium_dependent' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    //! Assembles the L2-projection operator
    virtual void assembleProjector() {
        std::cout << "Function "
                     "'assemble_projector' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    //! Releases all the memory allocated by this class
    void reset();

    //! Assembles the blocks of the linear system matrix
    void assembleSystem();

    //! Builds the linear system matrix from the blocks
    //! as an MFEM operator
    void buildSystemOp();

    //! Builds the linear system matrix from the blocks
    //! as a monolithic matrix
    void buildSystemMatrix();

    //! Builds the upper-triangle linear system matrix
    //! from the blocks as a monolithic matrix
    void buildSystemMatrixUpperTriangle();

    //! Assembles the target right-hand side
    void assembleRhs(mfem::BlockVector *) const;

    // for L2-projection
    //! Builds the L2-projection matrix
    void buildProjectorMatrix();

    //! Assembles the target right-hand side for L2-projection
    void assembleProjectionRhs(mfem::BlockVector *) const;

protected:

    virtual void assembleSourceXDiv(mfem::Vector&, double) const
    {
        std::cout << "Function "
                     "'assemble_source_xDiv' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    virtual void assembleXProjectionRhs(mfem::Vector&, mfem::Vector&,
                                          double) const {
        std::cout << "Function "
                     "'assemble_xProjection_rhs' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    //! Assembles the initial conditions
    void assembleICs(mfem::Vector&) const;

    void assembleICsXProjection(mfem::Vector&) const;

    void assembleSourceTGradXProjection(mfem::Vector&) const;
    void assembleSourceXProjection(mfem::Vector&, double) const;
    void assembleSourceTProjectionXDiv(mfem::Vector&) const;

    void addVector(const mfem::Array<int> &, const mfem::Vector&, mfem::Vector&) const;

    //! Applies boundary conditions on a matrix
    void applyBCs(mfem::SparseMatrix&) const;

    //! Applies boundary conditions on a vector
    void applyBCs(mfem::BlockVector&) const;

public:
    //! Returns the test case
    std::shared_ptr<heat::TestCases> getTestCase() const {
        return m_testCase;
    }
    
    //! Returns the temporal mesh
    mfem::Mesh* getTMesh() const {
        return m_tFespace->GetMesh();
    }
    
    //! Returns the spatial mesh
    mfem::Mesh* getXMesh() const {
        return m_xFespaces[0]->GetMesh();
    }

    //! Returns the temporal FE collection
    mfem::FiniteElementCollection* getTFec() const {
        return m_tFec;
    }

    //! Return the spatial FE collections
    mfem::Array<mfem::FiniteElementCollection*> getXFecs() const {
        return m_xFecs;
    }
    
    //! Returns the temporal FE space
    mfem::FiniteElementSpace* getTFespace() const {
        return m_tFespace;
    }
    
    //! Returns the spatial FE spaces
    mfem::Array<mfem::FiniteElementSpace*> getXFespaces() const {
        return m_xFespaces;
    }

    //! Returns the system operator
    mfem::BlockOperator* getHeatOp() const {
        return m_heatOp;
    }

    //! Returns the block offsets in the system operator,
    //! which has a block-structure
    mfem::Array<int> getBlockOffsets() const {
        return m_block_offsets;
    }

    //! Returns the monolithic system matrix
    mfem::SparseMatrix* getHeatMat() const {
        return m_heatMat;
    }

    //! Returns the projection matrix
    mfem::SparseMatrix* getHeatProjectionMat() const {
        return m_heatProjMat;
    }
    
    //! Prints some info
    void print() const {
        std::cout << "Degree in x and t: "
                  << m_deg << std::endl;
    }

    mfem::SparseMatrix* getSystemBlock11() const {
        return m_block00;
    }

    mfem::SparseMatrix* getSystemBlock12() const {
        return m_block01;
    }

    mfem::SparseMatrix* getSystemBlock21() const {
        return m_block10;
    }

    mfem::SparseMatrix* getSystemBlock22() const {
        return m_block11;
    }

protected:
    const nlohmann::json& m_config;
    
    int m_xndim, m_deg;
    std::shared_ptr<heat::TestCases> m_testCase;
    
    mfem::FiniteElementCollection* m_tFec = nullptr;
    mfem::FiniteElementSpace* m_tFespace = nullptr;
    
    mfem::Array<mfem::FiniteElementCollection*> m_xFecs;
    mfem::Array<mfem::FiniteElementSpace*> m_xFespaces;
    mfem::Array<int> m_block_offsets;
    
    mfem::SparseMatrix *m_tMass = nullptr;
    mfem::SparseMatrix *m_tStiff = nullptr;
    mfem::SparseMatrix *m_tGrad = nullptr;

    mfem::SparseMatrix *m_xMass1 = nullptr;
    mfem::SparseMatrix *m_xMass2 = nullptr;
    mfem::SparseMatrix *m_xStiff1 = nullptr;
    mfem::SparseMatrix *m_xStiff2 = nullptr;
    mfem::SparseMatrix *m_xGrad = nullptr;
    mfem::SparseMatrix *m_xDiv = nullptr;

    mfem::SparseMatrix *m_block00 = nullptr;
    mfem::SparseMatrix *m_block01 = nullptr;
    mfem::SparseMatrix *m_block10 = nullptr;
    mfem::SparseMatrix *m_block11 = nullptr;

    mfem::BlockOperator *m_heatOp = nullptr;
    mfem::SparseMatrix *m_heatMat = nullptr;
    mfem::SparseMatrix *m_heatProjMat = nullptr;
    
    mfem::Array<int> m_xEssBdrMarker;
    mfem::Array<int> m_xEssTdofList;
    mfem::Array<int> m_essTdofList;

    bool m_firstPass = true;
    mfem::SparseMatrix *m_block00MedIdp = nullptr;
    mfem::SparseMatrix *m_block01MedIdp = nullptr;
};

/**
 * @brief Least-squares space-time FEM for the heat equation,
 * with H1 discretisation in both time and space
 */
class LsqXtFemH1H1 : public LsqXtFEM
{
public:
    LsqXtFemH1H1 (const nlohmann::json& config,
                  std::shared_ptr<heat::TestCases>& testCase)
        : LsqXtFEM(config, testCase) {}

    void set(std::shared_ptr<mfem::Mesh>&,
             std::shared_ptr<mfem::Mesh>&) override;

    void assembleSystemMediumIndependent() override;
    void assembleSystemMediumDependent() override;

    void assembleProjector() override;

private:
    void assembleSourceXDiv(mfem::Vector&, double) const override;
    void assembleXProjectionRhs(mfem::Vector&,
                                mfem::Vector&, double) const override;
};

/**
 * @brief Least-squares space-time FEM for the heat equation,
 * with H1 discretisation in time and
 * Hdiv discretisation in space
 */class LsqXtFemH1Hdiv : public LsqXtFEM
{
public:
    LsqXtFemH1Hdiv (const nlohmann::json& config,
                        std::shared_ptr<heat::TestCases>& testCase)
        : LsqXtFEM(config, testCase) {}

    void set(std::shared_ptr<mfem::Mesh>&,
             std::shared_ptr<mfem::Mesh>&) override;

    void assembleSystemMediumIndependent() override;
    void assembleSystemMediumDependent() override;

    void assembleProjector() override;

private:
    void assembleSourceXDiv(mfem::Vector&, double) const override;
    void assembleXProjectionRhs(mfem::Vector&, mfem::Vector&,
                                  double) const override;
};

}

#endif // HEAT_DISCRETISATION_HPP
