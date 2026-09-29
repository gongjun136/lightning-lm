cmake_minimum_required(VERSION 3.16)
project(lightning_lm_docs NONE)

# Standalone documentation build: no ROS, compiler or application dependencies.
find_package(Doxygen REQUIRED COMPONENTS dot)
get_filename_component(LIGHTNING_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
configure_file(Doxyfile.in "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" @ONLY)
add_custom_target(docs ALL
        COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
        COMMAND "${CMAKE_COMMAND}"
                "-DWARN_LOG=${CMAKE_CURRENT_BINARY_DIR}/doxygen-warnings.log"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/check_warnings.cmake"
        WORKING_DIRECTORY "${LIGHTNING_SOURCE_DIR}"
        COMMENT "Generating Lightning-LM developer and algorithm manual"
        VERBATIM)
