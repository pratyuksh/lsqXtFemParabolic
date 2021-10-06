#ifndef PARDISO_HPP
#define PARDISO_HPP

#include "../core/config.hpp"

// PARDISO prototype.
#ifdef LIB_PARDISO
extern "C" void pardisoinit (void *, int    *, int *,
                             int  *, double *, int *);
extern "C" void pardiso     (void *, int    *, int *, int    *, int    *,
                             int  *, double *, int *, int    *, int    *,
                             int  *, int    *, int *, double *, double *,
                             int  *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *,
                                    int *, int *, int    *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_residual (int *, int * ,
                                  double *, int *, int *,
                                  double *, double *, double *,
                                  double *, double *);
#else // MKL Pardiso
extern "C" void PARDISO     (void *, int    *, int *, int    *, int    *,
                             int  *, double *, int *, int    *, int    *,
                             int  *, int    *, int *, double *, double *,
                             int  *);
#endif


/**
 * @brief Provides the function wrappers to use PARDISO solver,
 * available directly or from intel MKL library.
 */
class PardisoSolver
{
public:
    /**
     * @brief PardisoSolver Constructor
     * @param mtype Type of matrix, see PARDISO documentation
     */
    PardisoSolver (int mtype);

    /**
     * @brief PardisoSolver Constructor
     * @param config JSON config object
     */
    PardisoSolver (const nlohmann::json& config)
    {
        m_verbose = config["pardiso_verbose"];

        m_mtype = 11; // real unsymemtric matrix
        m_solver = 0; // sparse direct solver
    }

    //! Default destructor
    ~PardisoSolver ();

    /**
     * @brief Initializes the PARDISO solver
     * for a matrix given in the CSR format.
     * @param sizeA Dimension of the square-matrix, number of rows.
     * @param rowPtrA Row pointer in the CSR format
     * @param colIdA Column indices in the CSR format
     * @param dataA Matrix entries
     */
    void initialize (int sizeA, int *rowPtrA, int *colIdA, double *dataA);

    //! Calls the finalize routines in the PARDISO solver.
    void finalize();

    //! Factorize the linear system.
    void factorize ();

    /**
     * @brief Solves the linear system
     * @param b right-hand side vector
     * @param x solution vector
     */
    void solve (double *b, double *x);

private:
    int m_nnodes, m_nprocs;

    int m_mtype;
    int m_sizeA, m_nnz;
    int *m_rowPtrA, *m_colIdA; // CSR storage format
    double *m_dataA;

    int m_nrhs = 1;
    int m_solver;
    int m_maxfct, m_mnum, m_phase, m_error, m_verbose;

    int m_idum;
    double m_ddum;

    void *m_pt[64];
    int m_iparm[64];
    double m_dparm[64];

    bool m_finalized = true;
};

#endif /// PARDISO_HPP
