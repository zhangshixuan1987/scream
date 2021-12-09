# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test fields
  CreateUnitTest(field "field_tests.cpp" scream_field)
endif()
