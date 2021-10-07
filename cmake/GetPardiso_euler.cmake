find_library(PARDISO NAMES libpardiso700-GNU820-X86-64.so
  PATHS
    /cluster/home/oschenk
)

list(APPEND PARDISO gfortran)
