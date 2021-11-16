#ifndef SPARSE_HEAT_TEMPORAL_ASSEMBLY_HPP
#define SPARSE_HEAT_TEMPORAL_ASSEMBLY_HPP

#include "mfem.hpp"
using namespace mfem;

#include <memory>


namespace sparseHeat {

class TemporalBlockMatrixAssembler
{
public:
    //! Constructor with end-time,
    //! and minimum and maximum temporal mesh levels passed as arguments
    TemporalBlockMatrixAssembler(double T, int minLevel, int maxLevel);

    //! Initializes the assembler
    void initialize();

    //! Assembles the different blocks of the block-matrix
    virtual void assemble() = 0;

    //! Evaluates characteristics of the different meshes
    void evalMeshCharacteristics();

    //! Evaluates sizes and offsets of the blocks in the block-matrix
    void evalBlockSizesAndOffsets();

    //! Allocates memory for the block matrix member variable
    void allocateBlockMatrix();

    //! Returns the block sizes
    Array<int> getBlockSizes() {
        return m_blockSizes;
    }

    //! Returns the block offsets
    Array<int> getBlockOffsets() {
        return m_blockOffsets;
    }

    //! Returns the block-matrix
    std::shared_ptr<BlockMatrix> getBlockMatrix() {
         return m_blockMatrix;
    }
     
protected:
    double m_endTime;
    int m_numLevels, m_minLevel, m_maxLevel;

    Array<int> m_numMeshElements;
    Array<double> m_meshSizes;

    Array<int> m_blockSizes;
    Array<int> m_blockOffsets;

    std::shared_ptr<BlockMatrix> m_blockMatrix;
};


/**
 * @brief Assembler for the temporal canonical matrix
 */
class TemporalCanonicalMatrixAssembler
        : public TemporalBlockMatrixAssembler
{
public:
    TemporalCanonicalMatrixAssembler(double T, int minLevel, int maxLevel)
        : TemporalBlockMatrixAssembler(T, minLevel, maxLevel) {}

    void assemble() override;
};

/**
 * @brief Assembler for the temporal mass matrix
 */
class TemporalMassMatrixAssembler
        : public TemporalBlockMatrixAssembler
{
public:
    TemporalMassMatrixAssembler(double T, int minLevel, int maxLevel)
        : TemporalBlockMatrixAssembler(T, minLevel, maxLevel) {}

    void assemble() override;
};

/**
 * @brief Assembler for the temporal stiffness matrix
 */
class TemporalStiffnessMatrixAssembler
        : public TemporalBlockMatrixAssembler
{
public:
    TemporalStiffnessMatrixAssembler(double T, int minLevel, int maxLevel)
        : TemporalBlockMatrixAssembler(T, minLevel, maxLevel) {}

    void assemble() override;
};

/**
 * @brief Assembler for the temporal gradient matrix
 */
class TemporalGradientMatrixAssembler
        : public TemporalBlockMatrixAssembler
{
public:
    TemporalGradientMatrixAssembler(double T, int minLevel, int maxLevel)
        : TemporalBlockMatrixAssembler(T, minLevel, maxLevel) {}

    void assemble() override;
};

}

#endif // SPARSE_HEAT_TEMPORAL_ASSEMBLY_HPP
