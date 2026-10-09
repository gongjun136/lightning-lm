#pragma once

#include <string>

#include <pcl/conversions.h>
#include <pcl/io/pcd_io.h>

namespace lightning::pcd_io {

// World/map clouds already contain transformed XYZ. CloudCompare applies the
// PCD VIEWPOINT on import, so export these clouds with an identity viewpoint.
// Keep in-memory sensor metadata and every point field unchanged. Local sensor
// submaps should continue using their own frame and the regular PCL writer.
template <typename PointT>
int SaveWorldCloudBinaryCompressed(const std::string& filename,
                                   const pcl::PointCloud<PointT>& cloud) {
    pcl::PCLPointCloud2 serialized;
    pcl::toPCLPointCloud2(cloud, serialized);
    pcl::PCDWriter writer;
    return writer.writeBinaryCompressed(filename, serialized,
                                        Eigen::Vector4f::Zero(),
                                        Eigen::Quaternionf::Identity());
}

}  // namespace lightning::pcd_io
