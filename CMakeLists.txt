cmake_minimum_required(VERSION 3.16)

# Export editor compile settings by default, including after clean builds.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
        "Export compile_commands.json for C++ code navigation")

project(lightning_lm)

# Default to C++17
if (NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif ()

# set to release by default
set(CMAKE_BUILD_TYPE Release)

# -------------------------------------------------------------
# 导入第三方库，同时做一些全局编译选项配置
add_definitions("-DPCL_NO_PRECOMPILE")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fopenmp -fPIC -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -fopenmp -fPIC -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -g -fopenmp -fPIC -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -g")


if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(cmake/packages.cmake)
add_definitions(-DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/") # 设置项目根目录路径

install(DIRECTORY config
        DESTINATION lib/${PROJECT_NAME}/
)
install(PROGRAMS scripts/diagnose_localization.sh scripts/diagnose_localization.py
        DESTINATION lib/${PROJECT_NAME})
install(FILES scripts/diagnose_contract.json DESTINATION lib/${PROJECT_NAME})

# Keep the offline analyzer synchronized with the actual Safety producer. This
# check also runs in builds with BUILD_TESTING disabled; updating the manifest
# requires reviewing the parser and regression cases, not just accepting a hash.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(diagnostic_check "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_diagnostic_contract.py")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${diagnostic_check}" --print-dependencies
        OUTPUT_VARIABLE diagnostic_sources OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE diagnostic_status)
if(NOT diagnostic_status EQUAL 0)
    message(FATAL_ERROR "Cannot read diagnostic contract dependencies")
endif()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        ${diagnostic_sources} "${diagnostic_check}"
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/diagnose_contract.json")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${diagnostic_check}" --print-id
        OUTPUT_VARIABLE LIGHTNING_DIAGNOSTIC_CONTRACT OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE diagnostic_status ERROR_VARIABLE diagnostic_error)
if(NOT diagnostic_status EQUAL 0)
    message(FATAL_ERROR "${diagnostic_error}")
endif()
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE LIGHTNING_DIAGNOSTIC_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT LIGHTNING_DIAGNOSTIC_REVISION)
    set(LIGHTNING_DIAGNOSTIC_REVISION unavailable)
endif()
execute_process(COMMAND git status --porcelain --untracked-files=normal WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE diagnostic_dirty RESULT_VARIABLE diagnostic_git_status ERROR_QUIET)
set(LIGHTNING_DIAGNOSTIC_DIRTY 0)
if(diagnostic_dirty OR NOT diagnostic_git_status EQUAL 0)
    set(LIGHTNING_DIAGNOSTIC_DIRTY 1)
endif()
execute_process(COMMAND git rev-parse --git-path HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE diagnostic_git_head OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
execute_process(COMMAND git symbolic-ref -q HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE diagnostic_git_ref OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(diagnostic_git_head)
    get_filename_component(diagnostic_git_head "${diagnostic_git_head}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${diagnostic_git_head}")
endif()
if(diagnostic_git_ref)
    execute_process(COMMAND git rev-parse --git-path "${diagnostic_git_ref}" WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            OUTPUT_VARIABLE diagnostic_ref_path OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    get_filename_component(diagnostic_ref_path "${diagnostic_ref_path}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${diagnostic_ref_path}")
endif()
string(TIMESTAMP LIGHTNING_DIAGNOSTIC_BUILT_AT "%Y-%m-%dT%H:%M:%SZ" UTC)
configure_file(cmake/diagnostic_build_info.h.in generated/diagnostic_build_info.h @ONLY)

# Apply architecture-specific options only to Lightning sources.
if (LIGHTNING_CPU_COMPILE_OPTIONS)
    add_compile_options(${LIGHTNING_CPU_COMPILE_OPTIONS})
endif ()

add_subdirectory(src)

# Register the installed shared libraries in this scope so ament_package()
# generates the relocatable runtime library-path hook for local_setup.*.
ament_export_libraries(${PROJECT_NAME}.libs miao.core miao.utils)

option(BUILD_DOCS "Build the Doxygen developer and algorithm manual" OFF)
if (BUILD_DOCS)
    find_package(Doxygen REQUIRED COMPONENTS dot)
    find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
    set(LIGHTNING_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
            "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
            @ONLY
    )
    add_custom_target(docs ALL
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/maintenance/check_links.py"
            COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
            COMMAND "${CMAKE_COMMAND}"
                    "-DWARN_LOG=${CMAKE_CURRENT_BINARY_DIR}/doxygen-warnings.log"
                    -P "${CMAKE_CURRENT_SOURCE_DIR}/docs/check_warnings.cmake"
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/maintenance/prepare_html.py"
                    --html "${DOXYGEN_OUTPUT_DIRECTORY}/html"
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/maintenance/check_links.py"
                    --html "${DOXYGEN_OUTPUT_DIRECTORY}/html"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMENT "Generating Lightning-LM developer and algorithm manual"
            VERBATIM
    )
endif ()

ament_package()
