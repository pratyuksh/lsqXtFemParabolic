#include <gtest/gtest.h>

#include "../src/core/config.hpp"
#include "../src/mymfem/utilities.hpp"

using namespace mfem;


/**
 * @brief Unit test for the function to compute upper-triangle
 */
TEST(MfemUtil, getUpperTriangle)
{
    // CSR sparse matrix
    int sizeA = 8;
    int rowPtrA[9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
    int colIdA[20] = { 0,    2,       5, 6,
                          1, 2,    4,
                             2,             7,
                                   3,       6,
                          1,
                             2,       5,    7,
                          1,             6,
                             2,          6, 7 };
    double dataA[20] = { 7.0,      1.0,           2.0, 7.0,
                             -4.0, 8.0,      2.0,
                                   1.0,                     5.0,
                                        7.0,           9.0,
                             -4.0,
                                   7.0,           3.0,      8.0,
                              1.0,                    11.0,
                                  -3.0,                2.0, 5.0 };

    // create a Sparse matrix
    SparseMatrix trueA(rowPtrA, colIdA, dataA, sizeA, sizeA,
                   false, false, true);

    // get upper triangle
    SparseMatrix& UA = getUpperTriangle (trueA);

    // diagA = UA + UA.transpose() - trueA
    SparseMatrix diagA (*Add(UA, *Transpose(UA)));
    diagA.Add(-1, trueA);

    Vector trueDiagA, errDiagA;
    trueA.GetDiag(trueDiagA);
    diagA.GetDiag(errDiagA);
    errDiagA -= trueDiagA;

    double TOL = 1E-8;
    ASSERT_LE(errDiagA.Normlinf(), TOL);
}
