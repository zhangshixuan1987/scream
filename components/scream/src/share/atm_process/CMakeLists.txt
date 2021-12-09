set(ATM_PROC_SRC
  atmosphere_process.cpp
  atmosphere_process_group.cpp
  atmosphere_process_dag.cpp
)

add_library(atm_process ${ATM_PROC_SRC})
target_include_directories(atm_process PUBLIC ${SCREAM_SRC_DIR} ${SCREAM_BINARY_DIR}/src)
target_link_libraries(atm_process PUBLIC scream_io scream_grid scream_field scream_share)

# We have some issues with RDC and repeated libraries in the link line.
# It's a known issue, and nvcc_wrapper has a flag for handling this.
if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  if(CMAKE_VERSION VERSION_LESS "3.13.0")
    target_link_librarires(atm_proc PUBLIC "--remove-duplicate-link-files")
  else()
    target_link_options(atm_proc PUBLIC "--remove-duplicate-link-files")
  endif()
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
