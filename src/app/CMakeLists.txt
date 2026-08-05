add_executable(run_frontend_offline
        run_frontend_offline.cc
)

target_link_libraries(run_frontend_offline
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(run_frontend_online
        run_frontend_online.cc
)

target_link_libraries(run_frontend_online
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(run_slam_offline
        run_slam_offline.cc
)

target_link_libraries(run_slam_offline
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(run_slam_online
        run_slam_online.cc
)

target_link_libraries(run_slam_online
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(run_loc_offline
        run_loc_offline.cc
)

target_link_libraries(run_loc_offline
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(run_loc_online
        run_loc_online.cc
)

target_link_libraries(run_loc_online
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(convert_voxel_slam_map
        convert_voxel_slam_map.cc
)

target_link_libraries(convert_voxel_slam_map
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(build_solid_database
        build_solid_database.cc
)

target_link_libraries(build_solid_database
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(test_ui
        test_ui.cc
)

target_link_libraries(test_ui
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

add_executable(test_ui_init_quit
        test_ui_init_quit.cc
)

target_link_libraries(test_ui_init_quit
        ${PROJECT_NAME}.libs
        ${third_party_libs}
)

install(TARGETS run_frontend_offline run_frontend_online run_slam_offline run_slam_online run_loc_offline run_loc_online
        convert_voxel_slam_map build_solid_database
        DESTINATION lib/${PROJECT_NAME})
