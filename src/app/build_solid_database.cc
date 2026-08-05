#include <gflags/gflags.h>
#include <glog/logging.h>

#include "core/localization/solid_relocalizer.h"

DEFINE_string(config, "", "Lightning-LM YAML configuration");
DEFINE_string(map_path, "", "Lightning-LM map package directory");

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_alsologtostderr = true;
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_config.empty() || FLAGS_map_path.empty()) {
        LOG(ERROR) << "--config and --map_path are required";
        return 2;
    }
    return lightning::loc::SolidRelocalizer::BuildDatabase(
               FLAGS_config, FLAGS_map_path)
               ? 0
               : 1;
}
