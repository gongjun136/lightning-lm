#include <gtsam/geometry/Pose3.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>

#include <Eigen/Geometry>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct State {
  double timestamp;
  gtsam::Pose3 pose;
  gtsam::Vector6 odometry_variance;
};

struct LoopEdge {
  size_t history;
  size_t current;
  gtsam::Pose3 relative_pose;
};

std::vector<State> load_states(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open state file: " + path);

  std::vector<State> states;
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream stream(line);
    std::vector<double> values;
    double value;
    while (stream >> value) values.push_back(value);
    if (values.empty()) continue;
    if (values.size() != 26) throw std::runtime_error("expected 26 fields in state row");

    Eigen::Quaterniond quaternion(values[7], values[4], values[5], values[6]);
    quaternion.normalize();
    gtsam::Vector6 variance;
    for (size_t index = 0; index < 6; ++index) variance(index) = values[20 + index];
    states.push_back({
        values[0],
        gtsam::Pose3(gtsam::Rot3(quaternion.toRotationMatrix()),
                     gtsam::Point3(values[1], values[2], values[3])),
        variance,
    });
  }
  if (states.size() < 2) throw std::runtime_error("state file has fewer than two poses");
  return states;
}

LoopEdge load_single_edge(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open edge file: " + path);

  std::string session1, session2;
  size_t history, current;
  double tx, ty, tz, qx, qy, qz, qw;
  if (!(input >> session1 >> session2 >> history >> current >> tx >> ty >> tz >> qx >> qy >> qz >> qw)) {
    throw std::runtime_error("edge file does not contain one valid edge");
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("edge file contains more than one edge");

  Eigen::Quaterniond quaternion(qw, qx, qy, qz);
  quaternion.normalize();
  return {
      history,
      current,
      gtsam::Pose3(gtsam::Rot3(quaternion.toRotationMatrix()), gtsam::Point3(tx, ty, tz)),
  };
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: apply_voxel_loop_pose_graph STATE_FILE EDGE_FILE OUTPUT_TUM\n";
    return 2;
  }

  try {
    const std::vector<State> states = load_states(argv[1]);
    const LoopEdge loop = load_single_edge(argv[2]);
    if (loop.history >= states.size() || loop.current >= states.size()) {
      throw std::runtime_error("loop edge index is outside the state trajectory");
    }

    gtsam::NonlinearFactorGraph graph;
    gtsam::Values initial;
    for (size_t index = 0; index < states.size(); ++index) initial.insert(index, states[index].pose);

    gtsam::Vector6 prior_variance;
    prior_variance.setConstant(1e-9);
    graph.add(gtsam::PriorFactor<gtsam::Pose3>(
        0, states.front().pose, gtsam::noiseModel::Diagonal::Variances(prior_variance)));

    for (size_t index = 1; index < states.size(); ++index) {
      const gtsam::Pose3 measurement = states[index - 1].pose.between(states[index].pose);
      graph.add(gtsam::BetweenFactor<gtsam::Pose3>(
          index - 1, index, measurement,
          gtsam::noiseModel::Diagonal::Variances(states[index - 1].odometry_variance)));
    }

    gtsam::Vector6 loop_variance;
    loop_variance.setConstant(1e-4);
    graph.add(gtsam::BetweenFactor<gtsam::Pose3>(
        loop.history, loop.current, loop.relative_pose,
        gtsam::noiseModel::Diagonal::Variances(loop_variance)));

    const double initial_error = graph.error(initial);
    gtsam::ISAM2Params parameters;
    parameters.relinearizeThreshold = 0.01;
    parameters.relinearizeSkip = 1;
    gtsam::ISAM2 optimizer(parameters);
    optimizer.update(graph, initial);
    for (int iteration = 0; iteration < 5; ++iteration) optimizer.update();
    const gtsam::Values result = optimizer.calculateEstimate();
    const double final_error = graph.error(result);

    std::ofstream output(argv[3]);
    if (!output) throw std::runtime_error("cannot open output file");
    output << std::fixed << std::setprecision(9);
    double squared_translation_change = 0.0;
    double maximum_translation_change = 0.0;
    for (size_t index = 0; index < states.size(); ++index) {
      const gtsam::Pose3 pose = result.at<gtsam::Pose3>(index);
      const Eigen::Quaterniond quaternion(pose.rotation().matrix());
      const gtsam::Point3 translation = pose.translation();
      const double change = (translation - states[index].pose.translation()).norm();
      squared_translation_change += change * change;
      maximum_translation_change = std::max(maximum_translation_change, change);
      output << states[index].timestamp << ' ' << translation.x() << ' ' << translation.y() << ' '
             << translation.z() << ' ' << quaternion.x() << ' ' << quaternion.y() << ' '
             << quaternion.z() << ' ' << quaternion.w() << '\n';
    }

    std::cout << std::setprecision(12)
              << "poses=" << states.size() << '\n'
              << "loop_history=" << loop.history << '\n'
              << "loop_current=" << loop.current << '\n'
              << "graph_error_before=" << initial_error << '\n'
              << "graph_error_after=" << final_error << '\n'
              << "translation_change_rmse_m="
              << std::sqrt(squared_translation_change / states.size()) << '\n'
              << "translation_change_max_m=" << maximum_translation_change << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
