# Project
project(atlas)

cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_STANDARD 14)

# compile commands: by default, compile in Release mode
# type 'cmake -DCMAKE_BUILD_TYPE=Debug ..' in terminal for compiling in Debug mode
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# add libraries
add_subdirectory(commonutilities)

# locate external libraries
find_package(ZLIB REQUIRED)
find_package(Armadillo REQUIRED)

set(LINK_LIBRARIES
        coreLib
        commonutilities
        ${ZLIB_LIBRARIES}
        ${ARMADILLO_LIBRARIES}
        )

include_directories(        
        bamtools
        commonutilities/core
        commonutilities/core/IntegrationTests        
        GLF
        PopulationTools
        Simulations
        Recalibration
        Vcf
        Tests       
        ${ARMADILLO_INCLUDE_DIRS}
)

#######################
# Create core library
#######################
# assemble list of cpp and h files in all source folders
file(GLOB_RECURSE BAMTOOLS bamtools/*)
file(GLOB CORE_SOURCES *.h *.cpp
					   GLF/*.h GLF/*.cpp
                       PopulationTools/*.h PopulationTools/*.cpp
                       Simulations/*.h Simulations/*.cpp
                       Recalibration/*.h Recalibration/*.cpp
                       Vcf/*.h Vcf/*.cpp
                       Tests/*.h Tests/*.cpp
        )
list(APPEND CORE_SOURCES ${BAMTOOLS})

#Core library
add_library(coreLib STATIC ${CORE_SOURCES})

#######################
# Executable
#######################

# add executable to be built
add_executable(${PROJECT_NAME}
        atlas.cpp
        )

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
        ${LINK_LIBRARIES}
        )

# add the binary tree to the search path for include files
target_include_directories(${PROJECT_NAME} PUBLIC)

#######################
# Tests
#######################
