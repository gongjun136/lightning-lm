cmake_minimum_required(VERSION 3.8)

project(diagnostic_monitor_interfaces)

find_package(ament_cmake REQUIRED)
find_package(rosidl_default_generators REQUIRED)

rosidl_generate_interfaces(${PROJECT_NAME}
        "msg/NodeHeartbeat.msg"
        "msg/TopicCommHeader.msg"
)

# The owning project can also build this package as an embedded subdirectory.
# Make the generated IDL discoverable to other embedded interface targets.
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(diagnostic_monitor_interfaces_FOUND TRUE PARENT_SCOPE)
    set(diagnostic_monitor_interfaces_DIR
            "${CMAKE_CURRENT_BINARY_DIR}/rosidl_adapter/${PROJECT_NAME}/cmake"
            PARENT_SCOPE)
    set(diagnostic_monitor_interfaces_IDL_FILES
            "msg/NodeHeartbeat.idl;msg/TopicCommHeader.idl"
            PARENT_SCOPE)
endif()

ament_export_dependencies(rosidl_default_runtime)
ament_package()
