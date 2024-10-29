cmake_minimum_required(VERSION 3.14)

# Library
add_library(htslib SHARED IMPORTED GLOBAL)

set(HTSLIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/htslib")
set(INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")

if(APPLE)
set(LIB ${INSTALL_DIR}/lib/libhts.dylib)
else()
set(LIB ${INSTALL_DIR}/lib/libhts.so)
endif()

#inspired by https://gist.github.com/cauliyang/c6d69c43a8a05266c6d05bc2c5324f45:
set(flags --disable-gcs --disable-s3 --disable-plugins)

find_package(ZLIB REQUIRED)
find_package(LibLZMA)
  if(LIBLZMA_FOUND)
    include_directories(SYSTEM ${LIBLZMA_INCLUDE_DIRS})
  else()
    list(APPEND flags --disable-lzma)
  endif()

  find_package(CURL)
  if(CURL_FOUND)
    include_directories(SYSTEM ${CURL_INCLUDE_DIRS})
  else()
    list(APPEND flags --disable-libcurl)
  endif()

  find_package(BZip2)
  if(BZIP2_FOUND)
    include_directories(SYSTEM ${BZIP2_INCLUDE_DIRS})
  else()
    list(APPEND flags --disable-bz2)
  endif()

find_program(MAKE NAMES gmake make mingw32-make REQUIRED)
find_program(AUTORECONF NAMES autoreconf)

include(ExternalProject)
ExternalProject_Add(htslib-ext
  URL https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
  UPDATE_COMMAND ""
  BUILD_IN_SOURCE 1
  PREFIX             ${INSTALL_DIR}
  SOURCE_DIR         ${HTSLIB_DIR}
  CONFIGURE_COMMAND  ${AUTORECONF} -i ${HTSLIB_DIR} COMMAND ${HTSLIB_DIR}/configure --prefix=${INSTALL_DIR} ${flags}
  BUILD_COMMAND      ${MAKE}
  INSTALL_COMMAND    ${MAKE} install
  BUILD_BYPRODUCTS ${LIB}
  )

add_dependencies(htslib htslib-ext)
set_target_properties(htslib PROPERTIES IMPORTED_LOCATION ${LIB})
target_include_directories(htslib INTERFACE "${HTSLIB_DIR}")
target_link_libraries(htslib INTERFACE "${LIB}")
