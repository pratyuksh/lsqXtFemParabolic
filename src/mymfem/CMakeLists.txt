target_sources(MyMfem
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/base_observer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/nested_hierarchy.cpp
  #PRIVATE ${CMAKE_CURRENT_LIST_DIR}/my_bilinearForm_integrators.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/my_bilinearForms.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
)
