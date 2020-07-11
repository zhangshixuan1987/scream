# Process the values of options provided in ScreamOptions to establish values for
# other SCREAM configurations.
# Note: NONE of the following settings is a CACHE entry. The user has NO CONTROL
#       over the following settings.
# Note: these settings will error out if some user choices are not allowed.
#       E.g, the user may request to enable FPE with a pack size larger than 1,
#       which is not allowed.

###################################
#         Helper functions        #
###################################

# Shortcut function, to print a variable
function (print_var var)
  message ("  ${var}: ${${var}}")
endfunction ()

# Mini function to check that scream is configured with meaningful pack sizes
function (check_pack_size master_pack_size pack_size name)
  math (EXPR PACK_MODULO "${master_pack_size} % ${pack_size}")
  if ((pack_size GREATER master_pack_size) OR (NOT PACK_MODULO EQUAL 0))
    message (FATAL_ERROR "Invalid '${name}' size of ${pack_size}. Needs to be <= ${master_pack_size} and be a factor of it")
  endif()
endfunction ()

####################################
#        Process Pack sizes        #
####################################

if (SCREAM_POSSIBLY_NO_PACK)
  set (SCREAM_POSSIBLY_NO_PACK_SIZE 1)
else ()
  set (SCREAM_POSSIBLY_NO_PACK_SIZE ${SCREAM_PACK_SIZE})
endif ()

if (SCREAM_PACK_CHECK_BOUNDS)
  message (WARNING "Warning! Packs bounds check feature not yet implemented. This option has no effect.")
endif()

# Checks on pack sizes relative to the master one:
check_pack_size(${SCREAM_PACK_SIZE} ${SCREAM_SMALL_PACK_SIZE} "small pack")
# This one is an internal check, as the user cannot set SCREAM_POSSIBLY_NO_PACK_SIZE now.
check_pack_size(${SCREAM_PACK_SIZE} ${SCREAM_POSSIBLY_NO_PACK_SIZE} "possibly no pack")

####################################
#         FPE-related checks       #
####################################

# TODO: should we simply "turn off" fpes?
if (SCREAM_FPE)
  if (NOT ${SCREAM_PACK_SIZE} EQUAL 1)
    message (FATAL_ERROR "Error! Cannot enable FPEs if the pack size is larger than 1.")
  endif()
  if (CUDA_BUILD)
    message (FATAL_ERROR "Error! Cannot enable FPEs in a CUDA build.")
  endif ()
endif()

####################################
#      Process dynamics dycore     #
####################################

string (TOUPPER "${SCREAM_DYNAMICS_DYCORE}" SCREAM_DYNAMICS_DYCORE)
if (CIME_BUILD)
  if ("${SCREAM_DYN_TARGET}" STREQUAL "preqx_kokkos")
    if (NOT "${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME" AND
        NOT "${SCREAM_DYNAMICS_DYCORE}" STREQUAL "")
      message (WARNING "Warning! User requested SCREAM_DYN_TARGET=preqx_kokkos, but SCREAM_DYNAMICS_DYCORE=${SCREAM_DYNAMICS_DYCORE}\n"
                       "         Overriding SCREAM_DYNAMICS_DYCORE with 'HOMME'")
    endif ()
    set (SCREAM_DYNAMICS_DYCORE "HOMME")
  endif()
endif()

if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  set (SCREAM_HAS_HOMME TRUE CACHE INTERNAL "Signal that Homme is compiled")
elseif("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "NONE")
  message (STATUS "SCREAM_DYNAMICS_DYCORE set to 'NONE'. Scream won't enable any test/code that needs dynamics.\n")
elseif (NOT "${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  message (FATAL_ERROR "The dynamics dycore '${SCREAM_DYNAMICS_DYCORE}' is currently not supported."
                       "If you wish to disable dynamics, do not specify anything for SCREAM_DYNAMICS_DYCORE, or set it to 'NONE'.")
endif()

string(TOUPPER "${SCREAM_DYNAMICS_DYCORE}" SCREAM_DYNAMICS_DYCORE)
if (NOT ${SCREAM_DOUBLE_PRECISION})
  # Homme cannot handle single precision, for now. This causes tests to fail.
  # Fixing this requires adding a config parameter to homme, to switch between
  # single and double. That must be done in the upstream repo (E3SM), before
  # we can support it here.
  # So, for now, if Homme is the requested dyn dycore AND single precision is
  # requested, we disable dynamics, printing a warning.
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    message(FATAL_ERROR "ERROR! Homme dycore cannot be used in a Single Precision build.\n"
                        "       Either turn Homme off, or select Double Precision.")
  endif()
endif()

####################################
#      Process testing options     #
####################################

if (NOT SCREAM_LIB_ONLY AND NOT SCREAM_BASELINES_ONLY)
  enable_testing()
  include(CTest)
endif ()

file(MAKE_DIRECTORY ${SCREAM_TEST_DATA_DIR})

########################################
#      Print scream cmake settings     #
########################################

message ("\n******************** SCREAM Settings **********************")

print_var(SCREAM_DOUBLE_PRECISION)
print_var(SCREAM_MIMIC_GPU)
print_var(SCREAM_FPE)
print_var(SCREAM_PACK_SIZE)
print_var(SCREAM_SMALL_PACK_SIZE)
print_var(SCREAM_POSSIBLY_NO_PACK_SIZE)
print_var(SCREAM_LINK_FLAGS)
print_var(SCREAM_FPMODEL)
print_var(SCREAM_MPIRUN_EXE)
print_var(SCREAM_LIB_ONLY)
print_var(SCREAM_DYNAMICS_DYCORE)

message ("***********************************************************\n")
