cmake_minimum_required(VERSION 3.16...3.27)
project(Hadr09)

# ----------------------------------------------------------------------------
# Find Geant4
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

# ----------------------------------------------------------------------------
# Include directories
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${Geant4_INCLUDE_DIR}
)

# ----------------------------------------------------------------------------
# Main executable
set(MAIN_EXECUTABLE Hadr09)
file(GLOB MAIN_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_executable(${MAIN_EXECUTABLE} Hadr09.cc ${MAIN_SOURCES})
target_link_libraries(${MAIN_EXECUTABLE} ${Geant4_LIBRARIES} dl)
install(TARGETS ${MAIN_EXECUTABLE} DESTINATION bin)

# ----------------------------------------------------------------------------
# Build analysis plugin libraries
file(GLOB ANALYSIS_SOURCES ${PROJECT_SOURCE_DIR}/analyses/*.cc)

foreach(src ${ANALYSIS_SOURCES})
  get_filename_component(aname ${src} NAME_WE)
  add_library(${aname} SHARED ${src})
  target_include_directories(${aname} PRIVATE ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(${aname} ${Geant4_LIBRARIES})
  set_target_properties(${aname} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins
    OUTPUT_NAME "${aname}"
  )
  install(TARGETS ${aname} DESTINATION lib)
endforeach()

# ----------------------------------------------------------------------------
# Optional: Copy runtime scripts
set(Hadr09_SCRIPTS)
foreach(_script ${Hadr09_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()
