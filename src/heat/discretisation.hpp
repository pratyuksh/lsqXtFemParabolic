#ifndef HEAT_DISCRETISATION_HPP
#define HEAT_DISCRETISATION_HPP

#include "mfem.hpp"
using namespace mfem;

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
    virtual void set(std::shared_ptr<Mesh>&,
                     std::shared_ptr<Mesh>&) {
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

    //! Builds the linear system matrix from the blocks as an MFEM operator
    void buildSystemOp();

    //! Builds the linear system matrix from the blocks as a monolithic matrix
    void buildSystemMatrix();

    //! Builds the upper-triangle linear system matrix from the blocks as a monolithic matrix
    void buildSystemMatrixUpperTriangle();

    //! Assembles the target right-hand side
    void assembleRhs(BlockVector *) const;

    // for L2-projection
    //! Builds the L2-projection matrix
    void buildProjectorMatrix();

    //! Assembles the target right-hand side for L2-projection
    void assembleProjectionRhs(BlockVector *) const;

protected:

    virtual void assembleSourceXDiv(Vector&, double) const
    {
        std::cout << "Function "
                     "'assemble_source_xDiv' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    virtual void assembleXProjectionRhs(Vector&, Vector&,
                                          double) const {
        std::cout << "Function "
                     "'assemble_xProjection_rhs' "
                     "not implemented "
                     "for the class HeatLsqXtFEM!"
                  << std::endl;
        abort();
    }

    //! Assembles the initial conditions
    void assembleICs(Vector&) const;

    void assembleICsXProjection(Vector&) const;

    void assembleSourceTGradXProjection(Vector&) const;
    void assembleSourceXProjection(Vector&, double) const;
    void assembleSourceTProjectionXDiv(Vector&) const;

    void addVector(const Array<int> &, const Vector&, Vector&) const;

    //! Applies boundary conditions on a matrix
    void applyBCs(SparseMatrix&) const;

    //! Applies boundary conditions on a vector
    void applyBCs(BlockVector&) const;

public:
    //! Returns the test case
    inline std::shared_ptr<heat::TestCases> getTestCase() const {
        return m_testCase;
    }
    
    //! Returns the temporal mesh
    inline Mesh* getTMesh() const {
        return m_tFespace->GetMesh();
    }
    
    //! Returns the spatial mesh
    inline Mesh* getXMesh() const {
        return m_xFespaces[0]->GetMesh();
    }

    //! Returns the temporal FE collection
    inline FiniteElementCollection* getTFec() const {
        return m_tFec;
    }

    //! Return the spatial FE collections
    inline Array<FiniteElementCollection*> getXFecs() const {
        return m_xFecs;
    }
    
    //! Returns the temporal FE space
    inline FiniteElementSpace* getTFespace() const {
        return m_tFespace;
    }
    
    //! Returns the spatial FE spaces
    inline Array<FiniteElementSpace*> getXFespaces() const {
        return m_xFespaces;
    }

    //! Returns the system operator
    inline BlockOperator* getHeatOp() const {
        return m_heatOp;
    }

    //! Returns the block offsets in the system operator,
    //! which has a block-structure
    inline Array<int> getBlockOffsets() const {
        return m_block_offsets;
    }

    //! Returns the monolithic system matrix
    inline SparseMatrix* getHeatMat() const {
        return m_heatMat;
    }

    //! Returns the projection matrix
    inline SparseMatrix* getHeatProjectionMat() const {
        return m_heatProjMat;
    }
    
    //! Prints some info
    inline void print() const {
        std::cout << "Degree in x and t: "
                  << m_deg << std::endl;
    }

protected:
    const nlohmann::json& m_config;
    
    int m_xndim, m_deg;
    std::shared_ptr<heat::TestCases> m_testCase;
    
    FiniteElementCollection* m_tFec = nullptr;
    FiniteElementSpace* m_tFespace = nullptr;
    
    Array<FiniteElementCollection*> m_xFecs;
    Array<FiniteElementSpace*> m_xFespaces;
    Array<int> m_block_offsets;
    
    SparseMatrix *m_tMass = nullptr;
    SparseMatrix *m_tStiff = nullptr;
    SparseMatrix *m_tGrad = nullptr;

    SparseMatrix *m_xMass1 = nullptr;
    SparseMatrix *m_xMass2 = nullptr;
    SparseMatrix *m_xStiff1 = nullptr;
    SparseMatrix *m_xStiff2 = nullptr;
    SparseMatrix *m_xGrad = nullptr;
    SparseMatrix *m_xDiv = nullptr;

    SparseMatrix *m_block00 = nullptr;
    SparseMatrix *m_block01 = nullptr;
    SparseMatrix *m_block10 = nullptr;
    SparseMatrix *m_block11 = nullptr;

    BlockOperator *m_heatOp = nullptr;
    SparseMatrix *m_heatMat = nullptr;
    SparseMatrix *m_heatProjMat = nullptr;
    
    Array<int> m_xEssBdrMarker;
    Array<int> m_xEssTdofList;
    Array<int> m_essTdofList;

    bool m_firstPass = true;
    SparseMatrix *m_block00MedIdp = nullptr;
    SparseMatrix *m_block01MedIdp = nullptr;
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

    void set(std::shared_ptr<Mesh>&,
             std::shared_ptr<Mesh>&) override;

    void assembleSystemMediumIndependent() override;
    void assembleSystemMediumDependent() override;

    void assembleProjector() override;

private:
    void assembleSourceXDiv(Vector&, double) const override;
    void assembleXProjectionRhs(Vector&, Vector&, double) const override;
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

    void set(std::shared_ptr<Mesh>&,
             std::shared_ptr<Mesh>&) override;

    void assembleSystemMediumIndependent() override;
    void assembleSystemMediumDependent() override;

    void assembleProjector() override;

private:
    void assembleSourceXDiv(Vector&, double) const override;
    void assembleXProjectionRhs(Vector&, Vector&,
                                  double) const override;
};

}

#endif // HEAT_DISCRETISATION_HPP
