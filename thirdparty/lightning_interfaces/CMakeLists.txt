cmake_minimum_required(VERSION 3.16)

project(lightning)

find_package(ament_cmake REQUIRED)
find_package(rosidl_default_generators REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
if (NOT diagnostic_monitor_interfaces_FOUND)
    find_package(diagnostic_monitor_interfaces REQUIRED)
endif ()

rosidl_generate_interfaces(${PROJECT_NAME}
        "${CMAKE_SOURCE_DIR}:msg/DebugMessage.msg"
        "${CMAKE_SOURCE_DIR}:msg/FaultStatus.msg"
        "${CMAKE_SOURCE_DIR}:msg/LocalizationStatus.msg"
        "${CMAKE_SOURCE_DIR}:msg/PipelineDiagnostics.msg"
        "${CMAKE_SOURCE_DIR}:msg/VehiclePose.msg"
        "${CMAKE_SOURCE_DIR}:srv/LocCmd.srv"
        "${CMAKE_SOURCE_DIR}:srv/SaveMap.srv"
        DEPENDENCIES std_msgs geometry_msgs sensor_msgs nav_msgs diagnostic_monitor_interfaces
)

# When the common diagnostic interfaces are embedded in lightning_lm, they are
# not imported CMake targets.  Recreate the ordering/link edges that
# rosidl_generate_interfaces normally gets from find_package().
if (DIAGNOSTIC_MONITOR_INTERFACES_EMBEDDED)
    foreach (_rosidl_suffix
            rosidl_generator_c
            rosidl_generator_cpp
            rosidl_generator_py
            rosidl_typesupport_c
            rosidl_typesupport_cpp
            rosidl_typesupport_fastrtps_c
            rosidl_typesupport_fastrtps_cpp
            rosidl_typesupport_introspection_c
            rosidl_typesupport_introspection_cpp)
        set(_consumer_target "lightning__${_rosidl_suffix}")
        set(_dependency_target "diagnostic_monitor_interfaces__${_rosidl_suffix}")
        if (TARGET ${_consumer_target} AND TARGET ${_dependency_target})
            add_dependencies(${_consumer_target} ${_dependency_target})
            get_target_property(_consumer_type ${_consumer_target} TYPE)
            if (_consumer_type STREQUAL "INTERFACE_LIBRARY")
                target_link_libraries(${_consumer_target} INTERFACE ${_dependency_target})
            elseif (_rosidl_suffix STREQUAL "rosidl_generator_py")
                # rosidl_generator_py uses the plain signature internally.
                target_link_libraries(${_consumer_target} ${_dependency_target})
            elseif (NOT _consumer_type STREQUAL "UTILITY")
                target_link_libraries(${_consumer_target} PRIVATE ${_dependency_target})
            endif ()
        endif ()
    endforeach ()
endif ()

file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/lightning" CONTENT "")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lightning"
        DESTINATION share/ament_index/resource_index/packages)
install(FILES package.xml
        DESTINATION share/${PROJECT_NAME})
