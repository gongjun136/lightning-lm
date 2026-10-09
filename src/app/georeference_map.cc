#include "core/maps/tiled_map.h"
#include "common/pcd_io.h"
#include <pcl/io/pcd_io.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <yaml-cpp/yaml.h>
#include <glog/logging.h>

// Never modifies an existing map; rebuild relocalization databases separately.
int main(int argc,char** argv) {
    google::InitGoogleLogging(argv[0]);
    if (argc!=5) {
        std::cerr << "Usage: georeference_map source.pcd accepted_transform.yaml source_primary.tum NEW_output_dir\n";
        return 2;
    }
    try {
        namespace fs=std::filesystem;
        const fs::path output=fs::absolute(argv[4]);
        if (fs::exists(output)) throw std::invalid_argument("output directory must not exist");
        auto manifest=YAML::LoadFile(argv[2]);
        if (!manifest["accepted"].as<bool>() || manifest["convention"].as<std::string>()!="target_from_source")
            throw std::invalid_argument("require an accepted target_from_source transform");
        auto n=manifest["target_from_source"];
        auto t=n["translation_m"].as<std::vector<double>>(), q=n["quaternion_xyzw"].as<std::vector<double>>();
        if (t.size()!=3 || q.size()!=4) throw std::invalid_argument("transform dimensions");
        lightning::Quatd quat(q[3],q[0],q[1],q[2]); lightning::Vec3d trans(t[0],t[1],t[2]);
        if (!quat.coeffs().allFinite() || std::abs(quat.norm()-1)>1e-6 || !trans.allFinite())
            throw std::invalid_argument("invalid rigid transform");
        const lightning::SE3 transform(quat,trans);
        std::ifstream trajectory(argv[3]); double stamp,x,y,z,qx,qy,qz,qw;
        std::string line;
        bool found=false;
        while (std::getline(trajectory,line)) {
            if (line.find_first_not_of(" \t\r")==std::string::npos || line[line.find_first_not_of(" \t")]=='#') continue;
            std::istringstream row(line);
            found=static_cast<bool>(row>>stamp>>x>>y>>z>>qx>>qy>>qz>>qw);
            break;
        }
        if (!found) throw std::invalid_argument("source trajectory must contain a TUM pose");
        lightning::Quatd start_q(qw,qx,qy,qz); lightning::Vec3d start_t(x,y,z);
        if (!start_q.coeffs().allFinite() || std::abs(start_q.norm()-1)>1e-6 || !start_t.allFinite())
            throw std::invalid_argument("invalid source start pose");
        lightning::CloudPtr cloud(new lightning::PointCloudType);
        if (pcl::io::loadPCDFile(argv[1],*cloud)<0 || cloud->empty()) throw std::runtime_error("cannot read nonempty PCD");
        for (auto& p:*cloud) {
            if (!p.getVector3fMap().allFinite()) throw std::invalid_argument("nonfinite map point");
            p.getVector3fMap()=(transform*p.getVector3fMap().cast<double>()).cast<float>();
        }
        fs::create_directories(output/"tiled");
        if (lightning::pcd_io::SaveWorldCloudBinaryCompressed((output/"map.pcd").string(),*cloud)<0)
            throw std::runtime_error("cannot write transformed PCD");
        lightning::TiledMap map;
        if (!map.ConvertFromFullPCD(cloud,transform*lightning::SE3(start_q,start_t),(output/"tiled").string()) ||
            !fs::exists(output/"tiled"/"index.txt")) throw std::runtime_error("cannot export tiled map");
        fs::copy_file(argv[2],output/"georeference.yaml");
        std::cout << "Created " << output << "; rebuild SOLiD/BTC databases in this frame.\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
