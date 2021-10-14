target_sources(MyMfem
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mybilinearform.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/base_observer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/mesh_hierarchy.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities.cpp
)
