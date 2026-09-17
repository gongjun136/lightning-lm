find_package(glog REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Pangolin REQUIRED)
find_package(OpenGL REQUIRED)
find_package(pcl_conversions REQUIRED)
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(std_srvs REQUIRED)
find_package(OpenCV REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(rosbag2_cpp REQUIRED)
find_package(diagnostic_monitor_interfaces QUIET)
find_package(lightning QUIET)

# OMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

include(CheckCXXCompilerFlag)

set(LIGHTNING_CPU_PROFILE "AUTO" CACHE STRING
        "CPU optimization profile: AUTO, PORTABLE, NATIVE, or ORIN")
set_property(CACHE LIGHTNING_CPU_PROFILE PROPERTY STRINGS AUTO PORTABLE NATIVE ORIN)
string(TOUPPER "${LIGHTNING_CPU_PROFILE}" LIGHTNING_CPU_PROFILE_NORMALIZED)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" LIGHTNING_SYSTEM_PROCESSOR_NORMALIZED)

set(LIGHTNING_CPU_COMPILE_OPTIONS "")

function(lightning_select_first_supported_cpu_flag OUTPUT_VAR)
    foreach (CPU_FLAG IN LISTS ARGN)
        string(MAKE_C_IDENTIFIER "${CPU_FLAG}" CPU_FLAG_ID)
        set(SUPPORT_VAR "LIGHTNING_COMPILER_SUPPORTS_${CPU_FLAG_ID}")
        check_cxx_compiler_flag("${CPU_FLAG}" ${SUPPORT_VAR})
        if (${SUPPORT_VAR})
            set(${OUTPUT_VAR} "${CPU_FLAG}" PARENT_SCOPE)
            return()
        endif ()
    endforeach ()
    set(${OUTPUT_VAR} "" PARENT_SCOPE)
endfunction()

set(LIGHTNING_CPU_PROFILE_EFFECTIVE "${LIGHTNING_CPU_PROFILE_NORMALIZED}")
if (LIGHTNING_CPU_PROFILE_NORMALIZED STREQUAL "AUTO")
    if (CMAKE_CROSSCOMPILING)
        if (LIGHTNING_SYSTEM_PROCESSOR_NORMALIZED MATCHES "^(aarch64|arm64)$")
            set(LIGHTNING_CPU_PROFILE_EFFECTIVE "ORIN")
        elseif (LIGHTNING_SYSTEM_PROCESSOR_NORMALIZED MATCHES "^(x86_64|amd64|i[3-6]86)$")
            set(LIGHTNING_CPU_PROFILE_EFFECTIVE "X86")
        else ()
            set(LIGHTNING_CPU_PROFILE_EFFECTIVE "PORTABLE")
            message(WARNING
                    "No automatic CPU optimization is defined for cross target "
                    "${CMAKE_SYSTEM_PROCESSOR}; using portable CPU settings")
        endif ()
    else ()
        set(LIGHTNING_CPU_PROFILE_EFFECTIVE "NATIVE")
    endif ()
endif ()

if (LIGHTNING_CPU_PROFILE_EFFECTIVE STREQUAL "PORTABLE")
    # Keep the compiler's architecture baseline for portable artifacts.
elseif (LIGHTNING_CPU_PROFILE_EFFECTIVE STREQUAL "NATIVE")
    if (CMAKE_CROSSCOMPILING)
        message(FATAL_ERROR
                "LIGHTNING_CPU_PROFILE=NATIVE is invalid during cross compilation")
    endif ()

    if (LIGHTNING_SYSTEM_PROCESSOR_NORMALIZED MATCHES "^(aarch64|arm64)$")
        lightning_select_first_supported_cpu_flag(NATIVE_CPU_FLAG
                -mcpu=native -march=native)
    else ()
        lightning_select_first_supported_cpu_flag(NATIVE_CPU_FLAG -march=native)
    endif ()
    if (NATIVE_CPU_FLAG)
        list(APPEND LIGHTNING_CPU_COMPILE_OPTIONS "${NATIVE_CPU_FLAG}")
    elseif (LIGHTNING_SYSTEM_PROCESSOR_NORMALIZED MATCHES "^(x86_64|amd64|i[3-6]86)$")
        lightning_select_first_supported_cpu_flag(X86_CPU_FLAG
                -march=x86-64-v2 -msse4.2)
        if (X86_CPU_FLAG)
            list(APPEND LIGHTNING_CPU_COMPILE_OPTIONS "${X86_CPU_FLAG}")
        endif ()
    endif ()
    if (NOT LIGHTNING_CPU_COMPILE_OPTIONS)
        message(WARNING
                "The compiler does not support native CPU optimization; "
                "using portable CPU settings")
    endif ()
elseif (LIGHTNING_CPU_PROFILE_EFFECTIVE STREQUAL "X86")
    lightning_select_first_supported_cpu_flag(X86_CPU_FLAG
            -march=x86-64-v2 -msse4.2)
    if (X86_CPU_FLAG)
        list(APPEND LIGHTNING_CPU_COMPILE_OPTIONS "${X86_CPU_FLAG}")
    else ()
        message(WARNING
                "No x86 CPU optimization flag is supported; "
                "using portable CPU settings")
    endif ()
elseif (LIGHTNING_CPU_PROFILE_EFFECTIVE STREQUAL "ORIN")
    if (NOT LIGHTNING_SYSTEM_PROCESSOR_NORMALIZED MATCHES "^(aarch64|arm64)$")
        message(FATAL_ERROR
                "LIGHTNING_CPU_PROFILE=ORIN requires an ARM64 target, but "
                "CMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}")
    endif ()

    # GEACX2 uses Cortex-A78AE (Armv8.2-A). Prefer exact CPU tuning, then
    # progressively fall back to options supported by older ARM64 compilers.
    lightning_select_first_supported_cpu_flag(ORIN_CPU_FLAG
            -mcpu=cortex-a78ae -mcpu=cortex-a78 -march=armv8.2-a)
    if (ORIN_CPU_FLAG)
        list(APPEND LIGHTNING_CPU_COMPILE_OPTIONS "${ORIN_CPU_FLAG}")
    else ()
        message(WARNING
                "No GEACX2 CPU optimization flag is supported; "
                "using portable CPU settings")
    endif ()
else ()
    message(FATAL_ERROR
            "Unknown LIGHTNING_CPU_PROFILE=${LIGHTNING_CPU_PROFILE}. "
            "Expected AUTO, PORTABLE, NATIVE, or ORIN")
endif ()

message(STATUS "Lightning CPU profile: ${LIGHTNING_CPU_PROFILE_NORMALIZED} "
        "(effective: ${LIGHTNING_CPU_PROFILE_EFFECTIVE})")
message(STATUS "Lightning target processor: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Lightning cross compiling: ${CMAKE_CROSSCOMPILING}")
message(STATUS "Lightning CPU compile options: ${LIGHTNING_CPU_COMPILE_OPTIONS}")

include_directories(
        ${OpenCV_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${GLOG_INCLUDE_DIRS}
        ${Pangolin_INCLUDE_DIRS}
        ${GLEW_INCLUDE_DIRS}
        ${tf2_INCLUDE_DIRS}
        ${pcl_conversions_INCLUDR_DIRS}
        ${rclcpp_INCLUDE_DIRS}
        ${rosbag2_cpp_INCLUDE_DIRS}
        ${nav_msgs_INCLUDE_DIRS}
)

include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/livox_ros_driver/rosidl_generator_cpp
        ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/geosun_msgs/rosidl_generator_cpp
)

include_directories(
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/thirdparty
)


set(third_party_libs
        ${PCL_LIBRARIES}
        ${OpenCV_LIBS}
        ${Pangolin_LIBRARIES}
        glog::glog gflags
        ${yaml-cpp_LIBRARIES}
        ${pcl_conversions_LIBRARIES}
        tbb
        ${rosbag2_cpp_LIBRARIES}
)

