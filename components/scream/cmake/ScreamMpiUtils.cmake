# Detect the library that provides MPI
macro (GetMpiDistributionName DISTRO_NAME)
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} 
                  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryCompileOpenMPI.cpp
                  OUTPUT_QUIET ERROR_QUIET
                  RESULT_VARIABLE HAVE_OPENMPI_MPI)

  if (HAVE_OPENMPI_MPI)
    set (DISTRO_NAME "openmpi")
  else ()
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} 
                    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryCompileMPICH.cpp
                    OUTPUT_QUIET ERROR_QUIET
                    RESULT_VARIABLE HAVE_MPICH_MPI)

    if (HAVE_MPICH_MPI)
      set (DISTRO_NAME "mpich")
    else ()
      set (DISTRO_NAME "unknown")
    endif()
  endif()
endmacro ()

# Disable MPI cxx binding
macro (DisableMpiCxxBindings)
  # OpenMPI disable inclusion of mpicxx.h if the OMPI_SKIP_MPICXX macro
  # is defined, while MPICH checks the MPICH_SKIP_MPICXX macro.
  # To keep things tidy and avoid defining pointless macros, we first
  # check the MPI distribution, in order to decide what to define

  set (DISTRO_NAME)
  GetMpiDistributionName (DISTRO_NAME)

  if ("${DISTRO_NAME}" STREQUAL "openmpi")
    add_definitions (-DOMPI_SKIP_MPICXX)
  elseif ("${DISTRO_NAME}" STREQUAL "mpich")
    add_definitions (-DMPICH_SKIP_MPICXX)
  else ()
    message (FATAL_ERROR "Unsupported MPI distribution.")
  endif()
endmacro()

# Set the MPI cxx backend compiler var name
macro (SetMpiCxxBackendCompilerVarName)
  set (DISTRO_NAME)
  GetMpiDistributionName (DISTRO_NAME)

  if ("${DISTRO_NAME}" STREQUAL "openmpi")
    set (MPI_CXX_BACKEND_COMPILER_VAR_NAME  OMPI_CXX)
  elseif ("${DISTRO_NAME}" STREQUAL "mpich")
    set (MPI_CXX_BACKEND_COMPILER_VAR_NAME MPICH_CXX)
  else ()
    message (FATAL_ERROR "Unsupported MPI distribution.")
  endif()
endmacro()