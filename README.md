# Numerical schemes for Parabolic PDEs.
The primary target is to have a C++ serial implementation of the least-squares space-time scheme for parabolic PDEs, specifically the heat equation, with Pardiso linear system solver.

## Dependencies
1. MFEM-4.1, you can read more [info](https://mfem.org) and [download](https://mfem.org/download).
Follow the [instructions](https://mfem.org/building/) for the serial-build of MFEM.

1. Eigen library, [info](https://eigen.tuxfamily.org), download version - 3.3.7.

1. For writing the config files [JSON](https://github.com/nlohmann/json), download version - 3.7.0

1. The formatting library [fmt](https://fmt.dev/6.0.0).

1. Unit testing with [GoogleTest](https://github.com/google/googletest) framework.

Important directories:

- src: all the source and header files.
- tests: all unit tests are written here.
- config_files: JSON files specifying input parameters.
