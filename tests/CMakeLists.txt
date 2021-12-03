target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_pardiso.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh_point_locator.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_functionals.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_nested_hierarchy.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_my_bilinear_forms.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_spatial_assembly.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_spatial_assembly_H1Hdiv.cpp  
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_spatial_assembly_H1H1.cpp  
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_temporal_assembly.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_solution_handler.cpp  
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_discretisation.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_solver.cpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_heat_error_evaluator.cpp
)
