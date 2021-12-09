# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test grids
  CreateUnitTest(grid "grid_tests.cpp" scream_grid
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
