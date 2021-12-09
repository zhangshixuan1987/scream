set(FIELD_SRC
  field_alloc_prop.cpp
  field_identifier.cpp
  field_header.cpp
  field_layout.cpp
  field_tracking.cpp
  field.cpp
  field_group.cpp
  field_manager.cpp
  field_property_checks/field_nan_check.cpp
  field_property_checks/field_within_interval_check.cpp
)

add_library(scream_field ${FIELD_SRC})
target_link_libraries(scream_field PUBLIC scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
