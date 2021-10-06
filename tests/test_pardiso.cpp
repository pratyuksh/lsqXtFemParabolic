#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../src/pardiso/pardiso.hpp"


TEST(Pardiso, unsymmSolver)
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

    // create an Eigen Sparse matrix from rowPtrA, colIdA, dataA
    Eigen::SparseMatrix<double> matA(sizeA, sizeA);
    matA.reserve(20);
    for (int i=0; i<sizeA; i++) {
        for (int j=rowPtrA[i]; j<rowPtrA[i+1]; j++) {
            matA.insert(i,colIdA[j]) = dataA[j];
        }
    }
    matA.makeCompressed();

    // rhs
    Eigen::VectorXd b(sizeA);
    //b = Eigen::VectorXd::Constant(sizeA, 1);
    b = Eigen::VectorXd::Random(sizeA);

    // solve with Eigen
    Eigen::VectorXd true_x(sizeA);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparselu(matA);
    true_x = sparselu.solve(b);

    // solve with Pardiso
    int mtype = 11;
    Eigen::VectorXd x(sizeA);
    PardisoSolver pardisoSolver(mtype);
    pardisoSolver.initialize(sizeA, rowPtrA, colIdA, dataA);
    pardisoSolver.factorize();
    pardisoSolver.solve(b.data(), x.data());
    pardisoSolver.finalize();

    double TOL = 1E-8;
    ASSERT_LE((true_x - x).norm(), TOL);
}


TEST(Pardiso, symmSolver)
{
    // CSR sparse matrix
    int sizeA = 8;
    int rowPtrA[9] = { 0, 4, 7, 9, 11, 14, 16, 17, 18 };
    int colIdA[18] = { 0,    2,       5, 6,
                          1, 2,    4,
                             2,             7,
                               3,       6,
                                  4, 5, 6,
                                     5,    7,
                                        6,
                                           7 };
    double dataA[18] = { 7.0,      1.0,           2.0, 7.0,
                             -4.0, 8.0,           2.0,
                                   1.0,                     5.0,
                                        7.0,           9.0,
                                             5.0, -1.0, 5.0,
                                                  0.0,      5.0,
                                                       11.0,
                                                             5.0 };

    // create an Eigen Sparse matrix from rowPtrA, colIdA, dataA
    Eigen::SparseMatrix<double> matA(sizeA, sizeA);
    matA.reserve(28);
    for (int i=0; i<sizeA; i++) {
        for (int j=rowPtrA[i]; j<rowPtrA[i+1]; j++) {
            matA.insert(i,colIdA[j]) = dataA[j];
            if (i != colIdA[j]) {
                matA.insert(colIdA[j], i) = dataA[j];
            }
        }
    }
    matA.makeCompressed();

    // rhs
    Eigen::VectorXd b(sizeA);
    for (int i=0; i<sizeA; i++) {
        b(i) = i+1;
    }

    // solve with Eigen
    Eigen::VectorXd true_x(sizeA);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparselu(matA);
    true_x = sparselu.solve(b);
    //std::cout << true_x << std::endl;

    // solve with Pardiso
    int mtype = -2;
    Eigen::VectorXd x(sizeA);
    PardisoSolver pardisoSolver(mtype);
    pardisoSolver.initialize(sizeA, rowPtrA, colIdA, dataA);
    pardisoSolver.factorize();
    pardisoSolver.solve(b.data(), x.data());
    pardisoSolver.finalize();

    double TOL = 1E-8;
    ASSERT_LE((true_x - x).norm(), TOL);
}
