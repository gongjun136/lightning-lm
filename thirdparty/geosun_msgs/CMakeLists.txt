cmake_minimum_required(VERSION 3.16)

project(geosun_msgs)

find_package(rosidl_default_generators REQUIRED)
find_package(std_msgs REQUIRED)

rosidl_generate_interfaces(${PROJECT_NAME}
        "msg/PosRes.msg"
        DEPENDENCIES std_msgs
)

install(FILES resource/geosun_msgs
        DESTINATION share/ament_index/resource_index/packages)
install(FILES package.xml
        DESTINATION share/${PROJECT_NAME})
