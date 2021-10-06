#ifndef HEAT_UTILITIES_HPP
#define HEAT_UTILITIES_HPP

#include "mfem.hpp"
using namespace mfem;


/**
 * @brief Builds the solution in a space-time element
 * @param uXtSol coefficients of the discrete space-time solution
 * @param tShape FE shape vector in temporal dimension at time t
 * @param tVdofs Indices of the dofs in temporal dimension
 * @param uSol coefficients of the discrete spatial solution at time t
 */
void buildXSolFG(const Vector& uXtSol, const Vector& tShape,
                 Array<int> tVdofs, Vector& uSol);


#endif // HEAT_UTILITIES_HPP
