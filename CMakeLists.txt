project(atlas LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)

# add libraries (each containing an own CMakeLists.txt)
add_subdirectory(commonutilities)
add_subdirectory(VCF)
add_subdirectory(bamtools)
add_subdirectory(BAM)
add_subdirectory(GenomeTasks)
add_subdirectory(GenotypeLikelihoods)
add_subdirectory(GLF)
add_subdirectory(PopulationTools)
add_subdirectory(Simulations)
add_subdirectory(TestUtilities)

# add the executable
add_executable(${PROJECT_NAME} atlas.cpp)

# locate external libraries
find_package(ZLIB REQUIRED)
find_package(Armadillo REQUIRED)

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC VCFLib PopulationToolsLib BAMLib GenomeTasksLib GenotypeLikelihoodsLib GLFLib SimulationsLib TestsLib commonutilities bamtoolsLib ${ZLIB_LIBRARIES} ${ARMADILLO_LIBRARIES})

# add the binary tree to the search path for include files
target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS}
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/commonutilities"
        "${PROJECT_SOURCE_DIR}/bamtools"
        "${PROJECT_SOURCE_DIR}/VCF"
        "${PROJECT_SOURCE_DIR}/GenotypeLikelihoods"
        "${PROJECT_SOURCE_DIR}/BAM"
        "${PROJECT_SOURCE_DIR}/GenomeTasks"
        "${PROJECT_SOURCE_DIR}/GLF"
        "${PROJECT_SOURCE_DIR}/PopulationTools"
        "${PROJECT_SOURCE_DIR}/Simulations"
        "${PROJECT_SOURCE_DIR}/TestUtilities"
        )

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

# compile commands
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


################## google tests ####################

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
            "${gmock_SOURCE_DIR}/include"
            )
endif()

# Now simply link against gtest or gtest_main as needed. Eg
set(TEST_SOURCES BAM/tests/TGenomePosition_Test.cpp BAM/tests/TBed_Test.cpp BAM/tests/TBamFile_Test.cpp)
add_executable(testAtlas ${TEST_SOURCES})
target_link_libraries(testAtlas PUBLIC gtest_main gmock_main VCFLib PopulationToolsLib BAMLib GenomeTasksLib GenotypeLikelihoodsLib GLFLib SimulationsLib TestsLib commonutilities bamtoolsLib ${ZLIB_LIBRARIES} ${ARMADILLO_LIBRARIES})


target_include_directories(testAtlas PUBLIC ${ARMADILLO_INCLUDE_DIRS}
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/commonutilities"
        "${PROJECT_SOURCE_DIR}/bamtools"
        "${PROJECT_SOURCE_DIR}/VCF"
        "${PROJECT_SOURCE_DIR}/GenotypeLikelihoods"
        "${PROJECT_SOURCE_DIR}/BAM"
        "${PROJECT_SOURCE_DIR}/GenomeTasks"
        "${PROJECT_SOURCE_DIR}/GLF"
        "${PROJECT_SOURCE_DIR}/PopulationTools"
        "${PROJECT_SOURCE_DIR}/Simulations"
        "${PROJECT_SOURCE_DIR}/TestUtilities"
        )

add_test(NAME example_test COMMAND testAtlas)