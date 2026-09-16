cmake_minimum_required(VERSION 3.16)

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
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2 -fopenmp -fPIC -g -ggdb")
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

if (NOT lightning_FOUND)
    message(STATUS "lightning interface package not found; building embedded interfaces")
    find_package(ament_cmake_python REQUIRED)
    add_subdirectory(thirdparty/lightning_interfaces)
    _ament_cmake_python_register_environment_hook()
    include_directories(
            ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/lightning_interfaces/rosidl_generator_cpp
    )
endif ()

install(DIRECTORY config
        DESTINATION lib/${PROJECT_NAME}/
)

add_subdirectory(thirdparty/livox_ros_driver)
add_subdirectory(thirdparty/geosun_msgs)

# livox_ros_driver2 is embedded as a subdirectory, so its generated messages
# are installed into this package's prefix. Register the embedded interface
# package as well so generic ROS 2 tools (notably ros2 bag play) can discover it.
install(FILES thirdparty/livox_ros_driver/resource/livox_ros_driver2
        DESTINATION share/ament_index/resource_index/packages)
install(FILES thirdparty/livox_ros_driver/package.xml
        DESTINATION share/livox_ros_driver2)

add_subdirectory(src)

ament_package()
