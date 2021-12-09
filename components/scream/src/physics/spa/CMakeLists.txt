add_library(spa atmosphere_prescribed_aerosol.cpp)
target_link_libraries(spa physics_share atm_process scream_share)

add_subdirectory(tests)
