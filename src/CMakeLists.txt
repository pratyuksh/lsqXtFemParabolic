add_subdirectory(core)
add_subdirectory(pardiso)
add_subdirectory(mymfem)
add_subdirectory(heat)
add_subdirectory(sparse_heat)

target_sources(Heat
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/heat.cpp
)
