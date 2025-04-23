cmake_minimum_required(VERSION 3.16...3.27)
project(Hadr09)

# ----------------------------------------------------------------------------
# Options
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
option(WITH_YODA "Build with YODA analysis support" OFF)

# Force the linker to keep YODA even if not used in that binary
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")
set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -Wl,--no-as-needed")

# ----------------------------------------------------------------------------
# Find Geant4
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})

# ----------------------------------------------------------------------------
# Optional: Find YODA using yoda-config
if(WITH_YODA)
  find_program(YODA_CONFIG_EXECUTABLE NAMES yoda-config PATHS ENV PATH)

  if(NOT YODA_CONFIG_EXECUTABLE)
    message(FATAL_ERROR "yoda-config not found. Please load a YODA environment or set PATH correctly.")
  endif()

  execute_process(
    COMMAND ${YODA_CONFIG_EXECUTABLE} --cppflags
    OUTPUT_VARIABLE YODA_CPPFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  separate_arguments(YODA_CPPFLAGS)

  execute_process(
    COMMAND ${YODA_CONFIG_EXECUTABLE} --ldflags
    OUTPUT_VARIABLE YODA_LDFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  separate_arguments(YODA_LDFLAGS)

  message(STATUS "YODA cppflags: ${YODA_CPPFLAGS}")
  message(STATUS "YODA ldflags: ${YODA_LDFLAGS}")
endif()

# ----------------------------------------------------------------------------
# Include directories
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${Geant4_INCLUDE_DIR}
)

if(WITH_YODA)
  include_directories(${YODA_CPPFLAGS})
endif()

# ----------------------------------------------------------------------------
# Main executable
set(MAIN_EXECUTABLE Hadr09)
file(GLOB MAIN_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_executable(${MAIN_EXECUTABLE} Hadr09.cc ${MAIN_SOURCES})

target_link_libraries(${MAIN_EXECUTABLE} ${Geant4_LIBRARIES} dl)

if(WITH_YODA)
  target_compile_options(${MAIN_EXECUTABLE} PRIVATE ${YODA_CPPFLAGS})
  target_link_libraries(${MAIN_EXECUTABLE} ${YODA_LDFLAGS} YODA)
endif()


install(TARGETS ${MAIN_EXECUTABLE} DESTINATION bin)

# ----------------------------------------------------------------------------
# Build analysis plugin libraries
file(GLOB ANALYSIS_SOURCES ${PROJECT_SOURCE_DIR}/analyses/*.cc)

foreach(src ${ANALYSIS_SOURCES})
  get_filename_component(aname ${src} NAME_WE)
  add_library(${aname} SHARED ${src})
  target_include_directories(${aname} PRIVATE ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(${aname} ${Geant4_LIBRARIES})

  if(WITH_YODA)
    target_compile_options(${aname} PRIVATE ${YODA_CPPFLAGS})
    target_link_libraries(${aname} ${YODA_LDFLAGS} YODA)
  endif()

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
