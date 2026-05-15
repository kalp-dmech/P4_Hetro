CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

FIND_PACKAGE(deal.II 9.5.0 QUIET
  HINTS
    ${deal.II_DIR}
    ${DEAL_II_DIR}
    $ENV{DEAL_II_DIR}
    /home/me24m086/local/dealii-install
    /home/me24m086/local/dealii-install/lib/cmake/deal.II
    ../
    ../../
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

SET(SOURCE_DIR "./source/")
INCLUDE_DIRECTORIES(include)

SET(OUTPUT_DIR "./output")
file(MAKE_DIRECTORY ${OUTPUT_DIR})

SET(PROJECT_RUNTIME_LIBRARY_PATHS
  /home/me24m086/local/gcc-11/lib64
  /home/me24m086/local/dealii-install/lib
  /home/me24m086/local/petsc-install/lib
  /home/me24m086/local/p4est/lib
  /home/me24m086/local/openmpi/lib
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(project)

ADD_LIBRARY(lib
        ${SOURCE_DIR}Constructor.cc
        ${SOURCE_DIR}variable_constructor.cc
        ${SOURCE_DIR}make_grid.cc
        ${SOURCE_DIR}setup_qph.cc
        ${SOURCE_DIR}initialize_random_porosity.cc
        ${SOURCE_DIR}setup_constraints_elastic.cc
        ${SOURCE_DIR}setup_system_elastic.cc
        ${SOURCE_DIR}assemble_system_elastic.cc
        ${SOURCE_DIR}solve_linear_system_elastic.cc
        ${SOURCE_DIR}setup_constraints_damage.cc
        ${SOURCE_DIR}setup_system_damage.cc
        ${SOURCE_DIR}assemble_system_damage.cc
        ${SOURCE_DIR}H_plus.cc
        ${SOURCE_DIR}solve_linear_system_damage.cc
        ${SOURCE_DIR}check_convergence.cc
        ${SOURCE_DIR}output_results.cc
        ${SOURCE_DIR}update_history_field.cc
        ${SOURCE_DIR}load_disp_calculation.cc
        ${SOURCE_DIR}compute_energies.cc
        ${SOURCE_DIR}refine_grid.cc
        ${SOURCE_DIR}run.cc
        ${SOURCE_DIR}sed_data_post_processor.cc
        ${SOURCE_DIR}stress_post_processor.cc
        ${SOURCE_DIR}strain_post_processor.cc
        ${SOURCE_DIR}vector_post_processor.cc
        )
DEAL_II_SETUP_TARGET(lib)

ADD_EXECUTABLE(excecutable_output Phase_Field.cc)
DEAL_II_SETUP_TARGET(excecutable_output)
TARGET_LINK_LIBRARIES(excecutable_output lib)
SET_TARGET_PROPERTIES(excecutable_output PROPERTIES
  BUILD_RPATH "${PROJECT_RUNTIME_LIBRARY_PATHS}"
  INSTALL_RPATH "${PROJECT_RUNTIME_LIBRARY_PATHS}"
  )

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(run
  COMMAND excecutable_output
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Run with ${COMMAND_BUILD_TYPE} configuration")
