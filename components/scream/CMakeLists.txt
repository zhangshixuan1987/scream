# Detect whether we're building in standalone mode or not
if (NOT DEFINED PROJECT_NAME)
  cmake_minimum_required(VERSION 3.3)
  cmake_policy(SET CMP0057 NEW)
  set(CIME_BUILD FALSE)
else()
  set(CIME_BUILD TRUE)
endif()

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  # The new behavior for cmp0074 makes cmake use (rather than ignore)
  # any <PackageName>_ROOT env/cmake variable previously set.
  cmake_policy(SET CMP0074 NEW)
endif()

# Add the ./cmake folder to cmake path. Also add EKAT's cmake folder
list(APPEND CMAKE_MODULE_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake
     ${CMAKE_CURRENT_SOURCE_DIR}/ekat/cmake
)
include(EkatMpiUtils)

# We should avoid cxx bindings in mpi; they are already deprecated,
# and can cause headaches at link time, cause they require -lmpi_cxx
# (for openpmi; -lmpicxx for mpich) flag.
DisableMpiCxxBindings()

# cxx compiler and standard must be set *BEFORE* project is called
if (KOKKOS_ENABLE_CUDA)
  set (KOKKOS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/kokkos)
  include (EkatSetNvccWrapper)
  EkatSetNvccWrapper()
endif()
set (CMAKE_CXX_STANDARD 11)

#################################
#     Start of CMake Project    #
#################################

# 0) Declare the project
project(SCREAM CXX Fortran)

if (NOT CIME_BUILD)
  # Print the sha of the last commit (useful to double check which version was tested on CDash)
  execute_process (COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE LAST_GIT_COMMIT_SHA
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(LAST_GIT_COMMIT_SHA ${LAST_GIT_COMMIT_SHA} CACHE STRING "The sha of the last git commit.")
  message(STATUS "The sha of the last commit is ${LAST_GIT_COMMIT_SHA}")
endif()

# 1) Process Kokkos. We need this to be the first thing, since kokkos arch can
#    impact some of our defaults
set (Kokkos_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/kokkos)
include(Kokkos)

# 2) List and expose all SCREAM config options
include (ScreamOptions)

# 3) Process all the options (check the combination is supported, set internal configs, ...)
#    to determing the final scream settings.
#    Note: this also prints important settings for the user to see
include (ScreamSettings)

# 4) Set the some important directories, to be used across subfolders
set (SCREAM_BASE_DIR    ${CMAKE_CURRENT_SOURCE_DIR})
set (SCREAM_SRC_DIR     ${CMAKE_CURRENT_SOURCE_DIR}/src)
set (SCREAM_DATA_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/data)
set (SCREAM_F90_MODULES ${CMAKE_CURRENT_BINARY_DIR}/modules)

# 5) Proecess the EKAT package
set (YAMLCPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/yaml-cpp)
set (EKAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ekat)
include(Ekat)
if (SCREAM_DOUBLE_PRECISION)
  buildEkat("DOUBLE" PREFIX "SCREAM")
else()
  buildEkat("SINGLE" PREFIX "SCREAM")
endif()

# 6) Set compiler-specific flags
include(EkatSetCompilerFlags)
EkatSetCompilerFlags()

if (SCREAM_DOUBLE_PRECISION)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64 ")
  else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  endif()
endif()

# 7) Generate config files.
#    Note: fhis CPP macro must be done using add_definitions because it is used to determine
#    whether to include scream_config.h.
# TODO: do we really need this macro? SCREAM only works through cmake.
add_definitions(-DSCREAM_CONFIG_IS_CMAKE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/scream_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/scream_config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/scream_config.f.in ${CMAKE_CURRENT_BINARY_DIR}/src/scream_config.f)

# 8) Init scream TPLs info to whatever EKAT set.
set (SCREAM_TPL_LIBRARIES ${EKAT_TPL_LIBRARIES} CACHE INTERNAL "TPL libs that SCREAM targets can link to")
set (SCREAM_LINK_FLAGS ${KOKKOS_LDFLAGS_STR})
set (SCREAM_INCLUDE_DIRS
  ${SCREAM_SRC_DIR}
  ${EKAT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}/src
)

# 9) Declare baselines targets
if (NOT SCREAM_LIB_ONLY)
  add_custom_target(baseline)
  add_custom_target(baseline_cxx)
endif()

# 10) Process src and, possibly, test subfolders
add_subdirectory(src)
if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
