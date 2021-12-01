#include "utilities.hpp"

using namespace mfem;

// tVdofs gives the indices for reading the space-time solution uXtSol
// for the space-time element required.
// Assumes that uSol size is preset.
void buildXSolFG(const Vector& uXtSol,
                 const Vector& tShape,
                 Array<int> tVdofs,
                 Vector& uSol)
{
    int tNdofs = tShape.Size();
    int xdimV = uSol.Size();

    uSol = 0.0;
    for (int j=0; j<tNdofs; j++){
        int shift = tVdofs[j]*xdimV;
        for (int k=0; k<xdimV; k++) {
            uSol(k) += tShape(j)*uXtSol(k + shift);
        }
    }
}

// End of file
