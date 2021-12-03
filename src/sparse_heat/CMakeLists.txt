target_sources(LibSparseHeat
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spatial_assembly.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spatial_assembly_H1Hdiv.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spatial_assembly_H1H1.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/temporal_assembly.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/solution_handler.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/discretisation.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/discretisation_H1Hdiv.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/discretisation_H1H1.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/observer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spatial_error_evaluator.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/solver.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
)
