set(GRID_SRC
  abstract_grid.cpp
  mesh_free_grids_manager.cpp
  se_grid.cpp
  point_grid.cpp
  user_provided_grids_manager.cpp
  remap/abstract_remapper.cpp
)

add_library(scream_grid ${GRID_SRC})
target_link_libraries(scream_grid PUBLIC scream_field scream_share)

# We have some issues with RDC and repeated libraries in the link line.
# It's a known issue, and nvcc_wrapper has a flag for handling this.
if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  if(CMAKE_VERSION VERSION_LESS "3.13.0")
    target_link_librarires(scream_grid PUBLIC "--remove-duplicate-link-files")
  else()
    target_link_options(scream_grid PUBLIC "--remove-duplicate-link-files")
  endif()
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
