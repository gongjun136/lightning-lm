add_library(${PROJECT_NAME}.libs SHARED
        common/nav_state.cc
        common/constant.h
        common/options.cc
        common/params.cc

        io/yaml_io.cc
        io/file_io.cc

        core/lio/eskf.cc
        core/lio/laser_mapping.cc
        core/lio/multi_lidar_fusion.cc
        core/lio/pointcloud_preprocess.cc
        core/backend/backend_pipeline.cc
        core/backend/btc_loop_detector.cc
        core/backend/hierarchical_bundle_adjustment.cc
        core/backend/voxel_bundle_adjustment.cc
        core/backend/third_party/voxel_slam_btc/BTC.cpp
        core/loop_closing/loop_closing.cc

        core/g2p5/g2p5_map.cc
        core/g2p5/g2p5.cc
        core/g2p5/g2p5_subgrid.cc

        core/maps/tiled_map.cc
        core/maps/tiled_map_chunk.cc

        core/localization/localization_result.cc
        core/localization/lidar_loc/lidar_loc.cc
        core/localization/lidar_loc/pclomp/ndt_omp.cpp
        core/localization/lidar_loc/pclomp/voxel_grid_covariance_omp.cpp

        core/localization/pose_graph/pgo.cc
        core/localization/pose_graph/pgo_impl.cc
        core/localization/pose_graph/pose_extrapolator.cc
        core/localization/localization.cpp

        utils/timer.cc
        utils/pointcloud_utils.cc

        ui/pangolin_window.cc
        ui/pangolin_window_impl.cc
        ui/ui_car.cc
        ui/ui_cloud.cc
        ui/ui_trajectory.cc

        core/system/slam.cc
        core/system/loc_system.cc

        wrapper/bag_io.cc
)

ament_target_dependencies(
        ${PROJECT_NAME}.libs
        "rclcpp"
        "std_msgs"
        "geometry_msgs"
        "sensor_msgs"
        "nav_msgs"
        "std_srvs"
        "message_filters"
        "pcl_conversions"
        "tf2_ros"
        "tf2"
)

add_dependencies(${PROJECT_NAME}.libs
        livox_interfaces2__rosidl_typesupport_cpp
        ${PROJECT_NAME}__rosidl_typesupport_cpp
)

include_directories(${PROJECT_SOURCE_DIR}/src/core/miao)
add_subdirectory(core/miao)

target_link_libraries(${PROJECT_NAME}.libs
        ${third_party_libs}
        ${PCL_LIBRARIES}
        miao.core
        miao.utils
        livox_interfaces2__rosidl_typesupport_cpp
        ${PROJECT_NAME}__rosidl_typesupport_cpp
)

install(TARGETS ${PROJECT_NAME}.libs miao.core miao.utils
        DESTINATION lib)

if(BUILD_TESTING)
    add_executable(eskf_update_velocity_gate_test test/eskf_update_velocity_gate_test.cc)
    target_link_libraries(eskf_update_velocity_gate_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_executable(nav_state_eskf_model_test test/nav_state_eskf_model_test.cc)
    target_link_libraries(nav_state_eskf_model_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_executable(multi_lidar_fusion_test test/multi_lidar_fusion_test.cc)
    target_link_libraries(multi_lidar_fusion_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_executable(rear_axle_pose_test test/rear_axle_pose_test.cc)
    target_link_libraries(rear_axle_pose_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_executable(pointcloud_preprocess_test test/pointcloud_preprocess_test.cc)
    target_link_libraries(pointcloud_preprocess_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_executable(voxel_bundle_adjustment_test test/voxel_bundle_adjustment_test.cc)
    target_link_libraries(voxel_bundle_adjustment_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_executable(btc_loop_detector_test test/btc_loop_detector_test.cc)
    target_link_libraries(btc_loop_detector_test ${PROJECT_NAME}.libs ${third_party_libs})

    add_test(NAME eskf_update_velocity_gate_test COMMAND $<TARGET_FILE:eskf_update_velocity_gate_test>)
    add_test(NAME nav_state_eskf_model_test COMMAND $<TARGET_FILE:nav_state_eskf_model_test>)
    add_test(NAME multi_lidar_fusion_test COMMAND $<TARGET_FILE:multi_lidar_fusion_test>)
    add_test(NAME rear_axle_pose_test COMMAND $<TARGET_FILE:rear_axle_pose_test>)
    add_test(NAME pointcloud_preprocess_test COMMAND $<TARGET_FILE:pointcloud_preprocess_test>)
    add_test(NAME voxel_bundle_adjustment_test COMMAND $<TARGET_FILE:voxel_bundle_adjustment_test>)
    add_test(NAME btc_loop_detector_test COMMAND $<TARGET_FILE:btc_loop_detector_test>)
endif()

add_subdirectory(app)
