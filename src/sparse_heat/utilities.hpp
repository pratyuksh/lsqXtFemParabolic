#ifndef SPARSE_HEAT_UTILITIES_HPP
#define SPARSE_HEAT_UTILITIES_HPP

#include "mfem.hpp"


mfem::Array<int> evalBlockOffsets (mfem::Array<int>& blockSizes);

mfem::Array<int> evalTemporalBlockOffsets (int minLevel, int maxLevel);

mfem::Array<int> evalTemporalBlockSizes (int minLevel, int maxLevel);

mfem::Array<int> evalHierarchicalBlockSizes1D (int minLevel, int maxLevel);

mfem::Array<int> evalSpaceTimeBlockSizes (mfem::Array<int>& spatialBlockSizes,
                                          mfem::Array<int>& temporalBlockSizes);


mfem::Array<double> evalTemporalMeshSizes (double endTime,
                                           int minTemporalLevel,
                                           int numLevels);

void evalTemporalElIds (int finestTemporalElId,
                        mfem::Array<int>& temporalElIds);

void evalTemporalElLeftEdgeCoords (const mfem::Array<double> &temporalMeshSizes,
                                   const mfem::Array<int> &temporalElIds,
                                   mfem::Array<double> &temporalElLeftEdgeCoords);


double evalLeftHalfOfHatBasis(double z, double z0, double h);

double evalRightHalfOfHatBasis(double z, double z0, double h);

double affineTransform(double zLeft, double zRight, double eta);

#endif
