set(HAERO_SRCS
  atmosphere_aerosol.hpp
  haero_inputs_initializer.cpp
)

if (NOT CUDA_BUILD)
endif()

add_library(haero ${HAERO_SRCS})
target_include_directories(haero PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(haero PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_include_directories(haero SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS})
target_link_libraries(haero physics_share scream_share ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
 add_subdirectory(tests)
endif()
