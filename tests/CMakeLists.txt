target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_pardiso.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh_point_locator.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_nested_hierarchy.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_my_bilinear_forms.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_assembly.cpp  
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_functionals.cpp
)
