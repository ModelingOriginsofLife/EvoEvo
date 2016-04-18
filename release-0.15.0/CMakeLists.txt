
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the minimum required CMake version                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cmake_minimum_required(VERSION 2.6)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the project name                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
project(EVOEVO)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the package name                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(PACKAGE "\"MISEEM\"")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the version number                                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(VERSION_MAJOR 0)
set(VERSION_MINOR 15)
set(VERSION_PATCH 0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Set other options                                                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#set(WITH_GRAPHICS_CONTEXT 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Add definitions to CXX flags                                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ADD_DEFINITIONS(-std=c++11)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Set DEBUG and RELEASE flags                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g -pg -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3 -Wall -Wextra -pedantic")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the modules path                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the executables output path                                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(EXECUTABLE_OUTPUT_PATH ../build/bin)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the executables                                                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(CREATE_EXECUTABLE create)
add_executable(${CREATE_EXECUTABLE} src/create.cpp)

set(BOOTSTRAP_EXECUTABLE bootstrap)
add_executable(${BOOTSTRAP_EXECUTABLE} src/bootstrap.cpp)

set(RUN_EXECUTABLE run)
add_executable(${RUN_EXECUTABLE} src/run.cpp)

set(GENERATE_FIGURES_EXECUTABLE generate_figures)
add_executable(${GENERATE_FIGURES_EXECUTABLE} src/generate_figures.cpp)

set(UNITARY_TESTS_EXECUTABLE unitary_tests)
add_executable(${UNITARY_TESTS_EXECUTABLE} src/unitary_tests.cpp)

set(INTEGRATED_TESTS_EXECUTABLE integrated_tests)
add_executable(${INTEGRATED_TESTS_EXECUTABLE} src/integrated_tests.cpp)

set(DEEPEST_TREE_EXECUTABLE get_deepest_tree)
add_executable(${DEEPEST_TREE_EXECUTABLE} src/get_deepest_tree.cpp)

set(TROPHIC_GROUPS_EXECUTABLE trophic_groups_history)
add_executable(${TROPHIC_GROUPS_EXECUTABLE} src/trophic_groups_history.cpp)

set(PARAMETERS_EXECUTABLE recover_parameters)
add_executable(${PARAMETERS_EXECUTABLE} src/recover_parameters.cpp)

### SL post treatments ###

set(SL_FREQ_DEP_EXECUTABLE SL_frequency_dependence)
add_executable(${SL_FREQ_DEP_EXECUTABLE} src/SL_post_treatments/SL_frequency_dependence.cpp)

set(SL_CREATE_COMPETITION_EXECUTABLE SL_create_competition_experiment)
add_executable(${SL_CREATE_COMPETITION_EXECUTABLE} src/SL_post_treatments/SL_create_competition_experiment.cpp)

set(SL_RUN_COMPETITION_EXECUTABLE SL_run_competition_experiment)
add_executable(${SL_RUN_COMPETITION_EXECUTABLE} src/SL_post_treatments/SL_run_competition_experiment.cpp)

set(SL_HEATMAP_EXECUTABLE SL_heatmap)
add_executable(${SL_HEATMAP_EXECUTABLE} src/SL_post_treatments/SL_heatmap.cpp)

set(SL_CROSS_FEEDING_EXECUTABLE SL_cross_feeding)
add_executable(${SL_CROSS_FEEDING_EXECUTABLE} src/SL_post_treatments/SL_cross_feeding.cpp)

set(SL_BLACK_QUEEN_EXECUTABLE SL_black_queen)
add_executable(${SL_BLACK_QUEEN_EXECUTABLE} src/SL_post_treatments/SL_black_queen.cpp)

set(SL_DIVERSITY_EXECUTABLE SL_phylogenetic_diversity)
add_executable(${SL_DIVERSITY_EXECUTABLE} src/SL_post_treatments/SL_phylogenetic_diversity.cpp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Find and link external libraries                                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIR})
  target_link_libraries(${CREATE_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${RUN_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${DEEPEST_TREE_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${TROPHIC_GROUPS_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${PARAMETERS_EXECUTABLE} ${GSL_LIBRARIES})

  target_link_libraries(${SL_FREQ_DEP_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${SL_CREATE_COMPETITION_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${SL_RUN_COMPETITION_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${SL_HEATMAP_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${SL_CROSS_FEEDING_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${SL_BLACK_QUEEN_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${SL_DIVERSITY_EXECUTABLE} ${GSL_LIBRARIES})
endif(GSL_FOUND)

find_package(TBB REQUIRED)
if(TBB_FOUND)
  include_directories(${TBB_INCLUDE_DIR})
  target_link_libraries(${CREATE_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${RUN_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${DEEPEST_TREE_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${TROPHIC_GROUPS_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${PARAMETERS_EXECUTABLE} ${TBB_LIBRARIES})

  target_link_libraries(${SL_FREQ_DEP_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${SL_CREATE_COMPETITION_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${SL_RUN_COMPETITION_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${SL_HEATMAP_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${SL_CROSS_FEEDING_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${SL_BLACK_QUEEN_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${SL_DIVERSITY_EXECUTABLE} ${TBB_LIBRARIES})
endif(TBB_FOUND)

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(${CREATE_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${RUN_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${DEEPEST_TREE_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${TROPHIC_GROUPS_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${PARAMETERS_EXECUTABLE} ${ZLIB_LIBRARIES})

  target_link_libraries(${SL_FREQ_DEP_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${SL_CREATE_COMPETITION_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${SL_RUN_COMPETITION_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${SL_HEATMAP_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${SL_CROSS_FEEDING_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${SL_BLACK_QUEEN_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${SL_DIVERSITY_EXECUTABLE} ${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)

if(WITH_GRAPHICS_CONTEXT)
  find_package(SFML 2 REQUIRED graphics window system)
  if(SFML_FOUND)
    include_directories(${SFML_INCLUDE_DIR})
    target_link_libraries(${CREATE_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${RUN_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${DEEPEST_TREE_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${TROPHIC_GROUPS_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${PARAMETERS_EXECUTABLE} ${SFML_LIBRARIES})

    target_link_libraries(${SL_FREQ_DEP_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${SL_CREATE_COMPETITION_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${SL_RUN_COMPETITION_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${SL_HEATMAP_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${SL_CROSS_FEEDING_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${SL_BLACK_QUEEN_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${SL_DIVERSITY_EXECUTABLE} ${SFML_LIBRARIES})
  endif(SFML_FOUND)
endif(WITH_GRAPHICS_CONTEXT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create and link EvoEvo library                                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if(WITH_GRAPHICS_CONTEXT)
  add_library(
    EvoEvo
    src/lib/Cell.cpp
    src/lib/Cell.h
    src/lib/Enums.h
    src/lib/Environment.cpp
    src/lib/Environment.h
    src/lib/Genome.cpp
    src/lib/Genome.h
    src/lib/GraphicDisplay.cpp
    src/lib/GraphicDisplay.h
    src/lib/InheritedProteins.cpp
    src/lib/InheritedProteins.h
    src/lib/IntegratedTests.cpp
    src/lib/IntegratedTests.h
    src/lib/Macros.h
    src/lib/MutationEvent.cpp
    src/lib/MutationEvent.h
    src/lib/MutationVector.cpp
    src/lib/MutationVector.h
    src/lib/Node.cpp
    src/lib/Node.h
    src/lib/ODE.cpp
    src/lib/ODE.h
    src/lib/Parameters.cpp
    src/lib/Parameters.h
    src/lib/Population.cpp
    src/lib/Population.h
    src/lib/Prng.cpp
    src/lib/Prng.h
    src/lib/ReplicationReport.cpp
    src/lib/ReplicationReport.h
    src/lib/Simulation.cpp
    src/lib/Simulation.h
    src/lib/SpeciesList.cpp
    src/lib/SpeciesList.h
    src/lib/Statistics.cpp
    src/lib/Statistics.h
    src/lib/Structs.h
    src/lib/Tree.cpp
    src/lib/Tree.h
    src/lib/TrophicGroup.cpp
    src/lib/TrophicGroup.h
    src/lib/TrophicNetwork.cpp
    src/lib/TrophicNetwork.h
    src/lib/UnitaryTests.cpp
    src/lib/UnitaryTests.h
  )
elseif(NOT WITH_GRAPHICS_CONTEXT)
  add_library(
    EvoEvo
    src/lib/Cell.cpp
    src/lib/Cell.h
    src/lib/Enums.h
    src/lib/Environment.cpp
    src/lib/Environment.h
    src/lib/Genome.cpp
    src/lib/Genome.h
    src/lib/InheritedProteins.cpp
    src/lib/InheritedProteins.h
    src/lib/IntegratedTests.cpp
    src/lib/IntegratedTests.h
    src/lib/Macros.h
    src/lib/MutationEvent.cpp
    src/lib/MutationEvent.h
    src/lib/MutationVector.cpp
    src/lib/MutationVector.h
    src/lib/Node.cpp
    src/lib/Node.h
    src/lib/ODE.cpp
    src/lib/ODE.h
    src/lib/Parameters.cpp
    src/lib/Parameters.h
    src/lib/Population.cpp
    src/lib/Population.h
    src/lib/Prng.cpp
    src/lib/Prng.h
    src/lib/ReplicationReport.cpp
    src/lib/ReplicationReport.h
    src/lib/Simulation.cpp
    src/lib/Simulation.h
    src/lib/SpeciesList.cpp
    src/lib/SpeciesList.h
    src/lib/Statistics.cpp
    src/lib/Statistics.h
    src/lib/Structs.h
    src/lib/Tree.cpp
    src/lib/Tree.h
    src/lib/TrophicGroup.cpp
    src/lib/TrophicGroup.h
    src/lib/TrophicNetwork.cpp
    src/lib/TrophicNetwork.h
    src/lib/UnitaryTests.cpp
    src/lib/UnitaryTests.h
  )
endif(WITH_GRAPHICS_CONTEXT)

target_link_libraries(EvoEvo z tbb gsl gslcblas)

target_link_libraries(${CREATE_EXECUTABLE} EvoEvo)
target_link_libraries(${BOOTSTRAP_EXECUTABLE} EvoEvo)
target_link_libraries(${RUN_EXECUTABLE} EvoEvo)
target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} EvoEvo)
target_link_libraries(${UNITARY_TESTS_EXECUTABLE} EvoEvo)
target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} EvoEvo)
target_link_libraries(${DEEPEST_TREE_EXECUTABLE} EvoEvo)
target_link_libraries(${TROPHIC_GROUPS_EXECUTABLE} EvoEvo)
target_link_libraries(${PARAMETERS_EXECUTABLE} EvoEvo)

target_link_libraries(${SL_FREQ_DEP_EXECUTABLE} EvoEvo)
target_link_libraries(${SL_CREATE_COMPETITION_EXECUTABLE} EvoEvo)
target_link_libraries(${SL_RUN_COMPETITION_EXECUTABLE} EvoEvo)
target_link_libraries(${SL_HEATMAP_EXECUTABLE} EvoEvo)
target_link_libraries(${SL_CROSS_FEEDING_EXECUTABLE} EvoEvo)
target_link_libraries(${SL_BLACK_QUEEN_EXECUTABLE} EvoEvo)
target_link_libraries(${SL_DIVERSITY_EXECUTABLE} EvoEvo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build the Config.h file                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/modules/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
)
include_directories("${PROJECT_BINARY_DIR}")

