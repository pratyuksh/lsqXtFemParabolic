#include "utilities.hpp"

using namespace mfem;

// temporalVdofs gives the indices
// for reading the space-time solution,
// for the desired space-time element.
// Assumes that spatialSolution size is preset.
void buildSolutionAtSpecifiedTime(const Vector& spaceTimeSolution,
                                  const Vector& temporalShape,
                                  Array<int> temporalVdofs,
                                  Vector& spatialSolution)
{
    int temporalNumDofs = temporalShape.Size();
    int spatialNumDofs = spatialSolution.Size();

    spatialSolution = 0.0;
    for (int j=0; j<temporalNumDofs; j++)
    {
        int shift = temporalVdofs[j]*spatialNumDofs;
        for (int k=0; k<spatialNumDofs; k++) {
            spatialSolution(k)
                    += temporalShape(j)*spaceTimeSolution(k + shift);
        }
    }
}

// End of file
