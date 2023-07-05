cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Library
add_library(hts STATIC IMPORTED GLOBAL)

set(HTSLIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/htslib")
set(INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/htslib")
set(LIB "${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.so")

find_program(MAKE NAMES gmake make mingw32-make REQUIRED)
find_program(AUTORECONF NAMES autoreconf)

include(ExternalProject)
ExternalProject_Add(htslib
  GIT_REPOSITORY     "https://github.com/samtools/htslib.git"
  GIT_TAG            "master"
  UPDATE_DISCONNECTED true
  CONFIGURE_HANDLED_BY_BUILD true
  PREFIX             ${INSTALL_DIR}
  SOURCE_DIR         ${HTSLIB_DIR}
  CONFIGURE_COMMAND  ${AUTORECONF} -i ${HTSLIB_DIR} COMMAND ${HTSLIB_DIR}/configure --prefix=${INSTALL_DIR}
  BUILD_COMMAND      ${MAKE}
  INSTALL_COMMAND    ${MAKE} install
  BUILD_BYPRODUCTS ${LIB}
  )

add_dependencies(hts htslib)
set_target_properties(hts PROPERTIES IMPORTED_LOCATION ${LIB})
target_include_directories(hts INTERFACE "${HTSLIB_DIR}")
target_link_libraries(hts INTERFACE "${LIB}")
