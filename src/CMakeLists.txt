add_subdirectory(core)
add_subdirectory(pardiso)
add_subdirectory(mymfem)
add_subdirectory(heat)
add_subdirectory(sparse_heat)

target_sources(heat
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/heat.cpp
)

target_sources(sparse_heat
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/sparse_heat.cpp
)
