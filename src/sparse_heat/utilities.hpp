#ifndef SPARSE_HEAT_UTILITIES_HPP
#define SPARSE_HEAT_UTILITIES_HPP

#include "mfem.hpp"


mfem::Array<int> evalBlockOffsets (mfem::Array<int>& blockSizes);

mfem::Array<int> evalTemporalBlockOffsets (int minLevel, int maxLevel);

mfem::Array<int> evalTemporalBlockSizes (int minLevel, int maxLevel);

mfem::Array<int> evalHierarchicalBlockSizes1D (int minLevel, int maxLevel);

mfem::Array<int> evalSpaceTimeBlockSizes (mfem::Array<int>& spatialBlockSizes,
                                          mfem::Array<int>& temporalBlockSizes);


#endif
