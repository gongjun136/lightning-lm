#include "ui/pangolin_window.h"

#include <glog/logging.h>

#include <chrono>
#include <thread>

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;
    FLAGS_stderrthreshold = google::WARNING;
    google::ParseCommandLineFlags(&argc, &argv, true);

    lightning::ui::PangolinWindow ui;
    if (!ui.Init()) {
        LOG(ERROR) << "failed to init PangolinWindow";
        return 1;
    }
    if (ui.ShouldQuit()) {
        LOG(ERROR) << "PangolinWindow should not be quit immediately after init";
        return 2;
    }

    auto cloud = std::make_shared<lightning::PointCloudType>();
    cloud->reserve(200);
    for (int i = 0; i < 200; ++i) {
        lightning::PointType pt;
        pt.x = static_cast<float>(i % 20) - 10.0f;
        pt.y = static_cast<float>(i / 20) - 5.0f;
        pt.z = static_cast<float>((i % 7) - 3) * 0.1f;
        pt.intensity = static_cast<float>(i % 255);
        cloud->push_back(pt);
    }

    lightning::NavState state;
    state.pose_is_ok_ = true;
    state.timestamp_ = 1.0;
    state.rot_ = lightning::SO3();
    state.pos_ = lightning::Vec3d::Zero();
    ui.UpdateNavState(state);
    ui.UpdateScan(cloud, lightning::SE3());
    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    ui.Quit();
    if (!ui.ShouldQuit()) {
        LOG(ERROR) << "PangolinWindow should be quit after Quit()";
        return 3;
    }
    return 0;
}
