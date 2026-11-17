// 以测试专用点类型注入线搜索返回值，执行真实外层优化器；不向生产代码添加测试钩子。
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "ndt_omp.h"
#include "voxel_grid_covariance_omp_impl.hpp"
#include "ndt_omp_impl.hpp"

namespace {
int calls = 0, fail_at = 1;
double failure_value = std::numeric_limits<double>::quiet_NaN();
bool corrupt_hessian = false;
using Point = pcl::PointXYZRGB;
using Cloud = pcl::PointCloud<Point>;
using Ndt = pclomp::NormalDistributionsTransform<Point, Point>;
void Require(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
}

namespace pclomp {
// RGB 仅用于此独立测试，不覆盖生产定位器的 PointXYZIT 实例。
template <>
double NormalDistributionsTransform<Point, Point>::computeStepLengthMT(
    const Eigen::Matrix<double, 6, 1>& x, Eigen::Matrix<double, 6, 1>& direction,
    double, double, double, double&, Eigen::Matrix<double, 6, 1>&,
    Eigen::Matrix<double, 6, 6>& hessian, PointCloudSource& cloud) {
  ++calls;
  if (calls == fail_at) {
    if (corrupt_hessian) hessian(0, 0) = std::numeric_limits<double>::quiet_NaN();
    return failure_value;
  }
  if (fail_at == 0) return 0.0;  // 下一次正常调用应清除上次的异常状态。
  const Eigen::Matrix<double, 6, 1> trial = x + direction * 0.02;
  final_transformation_ = (Eigen::Translation3f(trial.head<3>().cast<float>()) *
      Eigen::AngleAxisf(trial(3), Eigen::Vector3f::UnitX()) *
      Eigen::AngleAxisf(trial(4), Eigen::Vector3f::UnitY()) *
      Eigen::AngleAxisf(trial(5), Eigen::Vector3f::UnitZ())).matrix();
  pcl::transformPointCloud(*input_, cloud, final_transformation_);
  return 0.02;
}
}

int main() {
  try {
    Cloud::Ptr cloud(new Cloud);
    for (int x = -2; x <= 2; ++x) for (int y = -2; y <= 2; ++y) for (int z = -1; z <= 1; ++z)
      for (int a : {-1, 1}) for (int b : {-1, 1}) for (int c : {-1, 1}) {
        Point p{}; p.x = x + a * .125f; p.y = y + b * .1875f; p.z = z + c * .21875f;
        cloud->push_back(p);
      }
    Ndt ndt;
    ndt.setResolution(1.0f); ndt.setNumThreads(1); ndt.setOulierRatio(.45);
    ndt.setNeighborhoodSearchMethod(pclomp::DIRECT7);
    ndt.AddTarget(cloud); ndt.ComputeTargetGrids(); ndt.setInputSource(cloud);
    ndt.setMaximumIterations(20); ndt.setTransformationEpsilon(.01);
    Eigen::Matrix4f guess = Eigen::Matrix4f::Identity(); guess(0, 3) = .08f;
    Cloud output;
    for (int failure_iteration : {1, 20}) for (double value : {
        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
      for (bool corrupt : {false, true}) {
        fail_at = failure_iteration; failure_value = value; corrupt_hessian = corrupt; calls = 0;
        ndt.align(output, guess);
        Require(calls == failure_iteration, "数值异常后仍继续调用线搜索");
        Require(ndt.getFinalNumIteration() == failure_iteration, "未保留实际尝试的迭代次数");
        Require(ndt.hasNumericalFailure() && !ndt.hasConverged(), "数值异常被误认为收敛");
        Require(std::isfinite(ndt.getTransformationProbability()) && ndt.getTransformationProbability() > .4,
                "未保留异常发生前可用的实际得分");
        Require(ndt.getFinalTransformation().allFinite() && ndt.getLastIncrementalTransformation().allFinite(),
                "非有限步长污染了位姿矩阵");
        Cloud expected; pcl::transformPointCloud(*cloud, expected, ndt.getFinalTransformation());
        Require(output.size() == expected.size(), "输出点云长度变化");
        for (std::size_t i = 0; i < output.size(); ++i)
          Require((output[i].getVector3fMap() - expected[i].getVector3fMap()).norm() < 1e-6f,
                  "保留的候选与输出点云不一致");
        fail_at = 0; calls = 0; ndt.align(output, guess);
        Require(!ndt.hasNumericalFailure() && ndt.hasConverged(), "下一次配准沿用了上次的异常状态");
      }
    std::cout << "NDT 数值异常退出回归通过：首轮／末轮、NaN／Inf、有限／异常 Hessian、复用恢复\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
