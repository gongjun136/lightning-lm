// ndt_line_search_test.cc
// 独立回归测试：恢复历史 NDT 步长策略，正常区间使用截断后的初始牛顿步。
// 仅生成合成点云，不访问文件/网络；不使用 assert；任一检查失败返回非零并打印中文说明。
// 体素键语义（分析方核实）：key = round(point * inverse_leaf_size)，分辨率 1.0 时整数为体素中心，半整数为边界。
// 预期：用例2确认不再进行内层步长细化；其余检查保留数值有效性和实际候选的一致性。
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <pcl/common/transforms.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include "ndt_omp.h"

namespace {

using PointT = pcl::PointXYZ;
using CloudT = pcl::PointCloud<PointT>;
using NdtT = pclomp::NormalDistributionsTransform<PointT, PointT>;
using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat6 = Eigen::Matrix<double, 6, 6>;

constexpr float kResolution = 1.0f;

// 仅用于暴露 protected 接口，不改变任何算法行为
class TestNdt : public NdtT {
 public:
  using NdtT::computeDerivatives;
  using NdtT::computeStepLengthMT;
};

int g_failures = 0;

bool Check(bool ok, const std::string& what) {
  std::printf("%s %s\n", ok ? "[通过]" : "[失败]", what.c_str());
  if (!ok) ++g_failures;
  return ok;
}

bool AllFinite(const Eigen::Matrix4f& m) {
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 4; ++c)
      if (!std::isfinite(m(r, c))) return false;
  return true;
}

// 确定性对称点簇：中心为整数格点，8 个立方角 ±d 与 6 个轴点 ±a_x/±a_y/±a_z。
// 所有偏移为二进制有理数，float 下 c±d 精确可表示 → 簇严格对称，均值恰为中心，协方差对角、三轴不等、非退化。
// 簇内最大偏移 0.21875，距最近体素边界（半整数）≥ 0.28。
CloudT::Ptr MakeClusters() {
  CloudT::Ptr c(new CloudT);
  const float d = 0.125f;
  const float ax = 0.21875f, ay = 0.1875f, az = 0.15625f;
  auto add = [&](float x, float y, float z) {
    PointT p;
    p.x = x;
    p.y = y;
    p.z = z;
    c->push_back(p);
  };
  for (int i = -2; i <= 2; ++i)
    for (int j = -2; j <= 2; ++j)
      for (int k = -1; k <= 1; ++k) {
        const float cx = static_cast<float>(i), cy = static_cast<float>(j), cz = static_cast<float>(k);
        for (int sx = -1; sx <= 1; sx += 2)
          for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2) add(cx + sx * d, cy + sy * d, cz + sz * d);
        add(cx + ax, cy, cz);
        add(cx - ax, cy, cz);
        add(cx, cy + ay, cz);
        add(cx, cy - ay, cz);
        add(cx, cy, cz + az);
        add(cx, cy, cz - az);
      }
  c->width = static_cast<uint32_t>(c->size());
  c->height = 1;
  c->is_dense = true;
  return c;
}

// 与 computeStepLengthMT 内部构造方式一致（float 截断、X-Y-Z 欧拉顺序）
Eigen::Matrix4f PoseMatrix(const Vec6& x) {
  return (Eigen::Translation<float, 3>(static_cast<float>(x(0)), static_cast<float>(x(1)),
                                      static_cast<float>(x(2))) *
          Eigen::AngleAxis<float>(static_cast<float>(x(3)), Eigen::Vector3f::UnitX()) *
          Eigen::AngleAxis<float>(static_cast<float>(x(4)), Eigen::Vector3f::UnitY()) *
          Eigen::AngleAxis<float>(static_cast<float>(x(5)), Eigen::Vector3f::UnitZ()))
      .matrix();
}

CloudT Transformed(const CloudT& source, const Vec6& x) {
  CloudT trans;
  pcl::transformPointCloud(source, trans, PoseMatrix(x));
  return trans;
}

// 按变体实际语义计算 round 体素键（分辨率 1.0）
std::vector<Eigen::Vector3i> KeysOf(const CloudT& c) {
  std::vector<Eigen::Vector3i> keys;
  keys.reserve(c.size());
  const float inv = 1.0f / kResolution;
  for (const PointT& p : c) {
    const Eigen::Vector3f v(p.x, p.y, p.z);
    keys.push_back((v * inv).array().round().cast<int>().matrix());
  }
  return keys;
}

bool SameKeys(const std::vector<Eigen::Vector3i>& a, const std::vector<Eigen::Vector3i>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[i]) return false;
  return true;
}

// 在位姿 x 处独立重算 score / gradient / Hessian（x 按值传入以满足非 const 引用参数）
double Derivs(TestNdt& ndt, const CloudT& source, Vec6 x, Vec6& g, Mat6& H, bool hess) {
  CloudT trans = Transformed(source, x);
  return ndt.computeDerivatives(g, H, trans, x, hess);
}

void Setup(TestNdt& ndt, const CloudT::Ptr& target, const CloudT::Ptr& source) {
  ndt.setResolution(kResolution);
  ndt.setNumThreads(1);
  ndt.setNeighborhoodSearchMethod(pclomp::DIRECT7);
  ndt.setStepSize(0.1);
  ndt.setTransformationEpsilon(1e-4);
  ndt.setMaximumIterations(50);
  ndt.AddTarget(target);
  ndt.ComputeTargetGrids();
  ndt.setInputSource(source);
}

// 运行一次 align：除作为端到端用例外，也用于初始化与分辨率相关的内部高斯常数，
// 使后续直接调用 computeDerivatives / computeStepLengthMT 时的状态与 align 内部一致。
Eigen::Matrix4f RunAlign(TestNdt& ndt, const Vec6& guess) {
  CloudT out;
  ndt.align(out, PoseMatrix(guess));
  return ndt.getFinalTransformation();
}

}  // namespace

// Precondition: include the NDT header under test (declaring pclomp::NormalDistributionsTransform)
// before this code. Call ndt_line_search_regression::run() from the existing test runner; true == pass.
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <pcl/common/transforms.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace ndt_line_search_regression {

typedef pcl::PointXYZ Pt;
typedef pcl::PointCloud<Pt> Cloud;
typedef Eigen::Matrix<double, 6, 1> Vec6;
typedef Eigen::Matrix<double, 6, 6> Mat6;

// Test subclass exposing the protected line-search helpers.
class Probe : public pclomp::NormalDistributionsTransform<Pt, Pt> {
 public:
  double trial(double a_l, double f_l, double g_l, double a_u, double f_u, double g_u, double a_t, double f_t,
               double g_t) {
    return this->trialValueSelectionMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t, f_t, g_t);
  }
  static Eigen::Matrix4f pose(const Vec6 &p) {
    return (Eigen::Translation<float, 3>(static_cast<float>(p(0)), static_cast<float>(p(1)),
                                         static_cast<float>(p(2))) *
            Eigen::AngleAxis<float>(static_cast<float>(p(3)), Eigen::Vector3f::UnitX()) *
            Eigen::AngleAxis<float>(static_cast<float>(p(4)), Eigen::Vector3f::UnitY()) *
            Eigen::AngleAxis<float>(static_cast<float>(p(5)), Eigen::Vector3f::UnitZ()))
        .matrix();
  }
  double evaluate(const Vec6 &x, Vec6 &g, Mat6 &h, Cloud &cloud) {
    Vec6 p = x;
    pcl::transformPointCloud(*this->input_, cloud, pose(p));
    return this->computeDerivatives(g, h, cloud, p, true);
  }
  double search(const Vec6 &x, Vec6 &dir, double step_init, double step_max, double step_min, double &score,
                Vec6 &g, Mat6 &h, Cloud &cloud) {
    return this->computeStepLengthMT(x, dir, step_init, step_max, step_min, score, g, h, cloud);
  }
};

struct Checker {
  int failures;
  Checker() : failures(0) {}
  void expect(bool ok, const char *what) {
    if (!ok) {
      ++failures;
      std::cerr << "NDT line-search regression failed: " << what << std::endl;
    }
  }
};

// Deterministic synthetic scene: jittered floor plus two walls. stride subsamples for the source.
inline Cloud::Ptr makeCloud(int stride) {
  Cloud::Ptr c(new Cloud);
  unsigned int s = 12345u;
  int n = 0;
  for (int i = 0; i <= 80; ++i) {
    for (int j = 0; j <= 80; ++j) {
      const float u = -4.0f + 0.1f * static_cast<float>(i);
      const float v = -4.0f + 0.1f * static_cast<float>(j);
      const float w = 0.0375f * static_cast<float>(j);
      for (int k = 0; k < 3; ++k) {
        s = s * 1664525u + 1013904223u;
        const float jit = static_cast<float>((s >> 8) & 0xFFFFu) / 65535.0f * 0.04f - 0.02f;
        Pt p;
        if (k == 0) {
          p.x = u; p.y = v; p.z = jit;
        } else if (k == 1) {
          p.x = 4.0f + jit; p.y = u; p.z = w;
        } else {
          p.x = u; p.y = 4.0f + jit; p.z = w;
        }
        if (n++ % stride == 0) c->push_back(p);
      }
    }
  }
  return c;
}

// Runs one line search from the evaluated state (x, score0, g0, h0), checks the postconditions and
// returns alpha. Core invariant: returned alpha, score, gradient, Hessian and trans_cloud all describe
// the same actually evaluated pose x + dir * alpha.
inline double runAndCheck(Probe &ndt, Checker &c, const Vec6 &x, const Vec6 &dir_in, double step_init,
                          double step_max, double step_min, double score0, const Vec6 &g0, const Mat6 &h0,
                          double *score_out) {
  Vec6 dir = dir_in;  // may be reversed in place by the line search
  double score = score0;
  Vec6 g = g0;
  Mat6 h = h0;
  Cloud cloud;
  const double alpha = ndt.search(x, dir, step_init, step_max, step_min, score, g, h, cloud);
  if (score_out) *score_out = score;

  c.expect(std::isfinite(alpha), "alpha is finite");
  c.expect(alpha > 0.0, "alpha is a real positive step, not an artificial zero");
  c.expect(alpha >= step_min && alpha <= step_max, "alpha lies within [step_min, step_max]");
  c.expect(std::isfinite(score) && g.allFinite() && h.allFinite(), "score, gradient and Hessian are finite");
  c.expect(h.norm() > 0.0, "selected candidate has a nonzero Hessian");
  if (!std::isfinite(alpha)) return alpha;

  Vec6 g2;
  Mat6 h2;
  Cloud cloud2;
  const Vec6 x_ret = x + dir * alpha;
  const double score2 = ndt.evaluate(x_ret, g2, h2, cloud2);
  const double tol = 1e-6;
  c.expect(std::fabs(score - score2) <= tol * (1.0 + std::fabs(score2)), "score matches returned alpha");
  c.expect((g - g2).norm() <= tol * (1.0 + g2.norm()), "gradient matches returned alpha");
  c.expect((h - h2).norm() <= tol * (1.0 + h2.norm()), "Hessian matches returned alpha");
  bool same = cloud.size() == cloud2.size();
  for (std::size_t i = 0; same && i < cloud.size(); ++i) {
    same = std::fabs(cloud[i].x - cloud2[i].x) <= 1e-4f && std::fabs(cloud[i].y - cloud2[i].y) <= 1e-4f &&
           std::fabs(cloud[i].z - cloud2[i].z) <= 1e-4f;
  }
  c.expect(same, "trans_cloud matches returned alpha");
  return alpha;
}

inline bool run() {
  Checker c;

  // 1. Degenerate interpolation, synthetic endpoints (characterises the unguarded helper; these pass
  //    before and after the fix and show why the caller must validate the candidate).
  {
    Probe helper;
    // Identical endpoints after clamping to step_min (a_t == a_l) with f_t > f_l: 0/0 in Case 1.
    const double d1 = helper.trial(0.005, 1.0, -1.0, 0.0, 0.0, -1.0, 0.005, 2.0, -1.0);
    c.expect(!std::isfinite(d1), "identical-endpoint interpolation is nonfinite in the raw helper");
    // Case 4 with negative discriminant: z == 0, g_t * g_u > 0.
    const double d2 = helper.trial(0.0, 0.0, -1.0, 2.0, -2.5, -1.0, 1.0, -1.5, -2.0);
    c.expect(!std::isfinite(d2), "negative-discriminant interpolation is nonfinite in the raw helper");
  }

  // 2. Synthetic scene. One public align() call initialises all internal state.
  Probe ndt;
  Cloud::Ptr target = makeCloud(1);
  Cloud::Ptr source = makeCloud(7);
  ndt.setResolution(1.0f);
  ndt.setStepSize(0.1);
  ndt.setTransformationEpsilon(0.01);
  ndt.setMaximumIterations(1);
  ndt.setNumThreads(1);
  ndt.AddTarget(target);
  ndt.ComputeTargetGrids();
  ndt.setInputSource(source);
  Cloud aligned;
  ndt.align(aligned, Eigen::Matrix4f::Identity());

  Vec6 x0;
  x0 << 0.3, -0.2, 0.05, 0.0, 0.0, 0.02;
  Vec6 g0;
  Mat6 h0;
  Cloud cloud0;
  const double score0 = ndt.evaluate(x0, g0, h0, cloud0);
  c.expect(std::isfinite(score0) && g0.allFinite() && h0.allFinite() && g0.norm() > 0.0,
           "setup: finite, non-stationary starting state");
  if (c.failures) return false;
  const Vec6 gdir = g0.normalized();

  // 3. Legacy positive-interval policy: even an oversized trial retains its clamped step.
  runAndCheck(ndt, c, x0, gdir, 2.0, 2.0, 0.005, score0, g0, h0, 0);

  // 4. Clamp stagnation. With step_min == step_max every proposal clamps onto the evaluated step,
  //    so the search must return exactly that step with a consistent state.
  {
    const double a_small = runAndCheck(ndt, c, x0, gdir, 0.005, 0.005, 0.005, score0, g0, h0, 0);
    c.expect(a_small == 0.005, "stagnation at the lower bound returns exactly step_min");
    const double a_large = runAndCheck(ndt, c, x0, gdir, 1.5, 1.5, 1.5, score0, g0, h0, 0);
    c.expect(a_large == 1.5, "stagnation at an overshooting clamped step returns exactly that step");
    const double a_below = runAndCheck(ndt, c, x0, gdir, 1e-6, 2.0, 0.005, score0, g0, h0, 0);
    c.expect(a_below >= 0.005, "initial step below step_min is clamped up, never zero");
  }

  // 5. Sweep with repeated calls: every call from the same inputs must be finite, bounded,
  //    state-consistent and reproducible (no state leaking between calls).
  {
    const double inits[5] = {0.005, 0.05, 0.4, 1.0, 2.0};
    for (int d = 0; d < 8; ++d) {
      Vec6 dir;
      if (d == 0) dir = gdir;
      else if (d == 1) dir = -gdir;
      else dir = Vec6::Unit(d - 2);
      if (g0.dot(dir) == 0.0) continue;  // pre-existing exact-zero-slope early return, out of scope
      for (int k = 0; k < 5; ++k) {
        double s1 = 0.0, s2 = 0.0;
        const double a1 = runAndCheck(ndt, c, x0, dir, inits[k], 2.0, 0.005, score0, g0, h0, &s1);
        const double a2 = runAndCheck(ndt, c, x0, dir, inits[k], 2.0, 0.005, score0, g0, h0, &s2);
        c.expect(std::fabs(a1 - a2) <= 1e-9, "repeated call returns the same alpha");
        c.expect(std::fabs(s1 - s2) <= 1e-6 * (1.0 + std::fabs(s1)), "repeated call returns the same score");
      }
    }
  }

  // 6. Nonfinite input direction: must signal NaN (not zero, not a finite step) without evaluating
  //    or modifying outputs.
  {
    Vec6 bad = gdir;
    bad(2) = std::numeric_limits<double>::quiet_NaN();
    double score = score0;
    Vec6 g = g0;
    Mat6 h = h0;
    Cloud cloud;
    const double alpha = ndt.search(x0, bad, 0.4, 2.0, 0.005, score, g, h, cloud);
    c.expect(alpha != alpha, "nonfinite direction yields NaN alpha, not a zero or finite step");
    c.expect(score == score0 && g == g0 && h == h0 && cloud.empty(), "nonfinite direction leaves outputs untouched");
  }

  // 7. The object is still usable after the rejected call.
  runAndCheck(ndt, c, x0, gdir, 0.4, 2.0, 0.005, score0, g0, h0, 0);

  return c.failures == 0;
}

}  // namespace ndt_line_search_regression

int main() {
  if (!ndt_line_search_regression::run()) return 1;
  // 目标与源为同一确定性点集（分别生成两份，避免共享指针），单位变换为已知最优
  CloudT::Ptr target = MakeClusters();
  CloudT::Ptr source = MakeClusters();
  TestNdt ndt;
  Setup(ndt, target, source);

  // 用例 1：端到端，从小偏移初值收敛到单位变换（健全性，不区分修复前后）
  {
    Vec6 guess;
    guess << 0.08, -0.06, 0.05, 0.010, -0.008, 0.012;
    const Eigen::Matrix4f T = RunAlign(ndt, guess);
    const Eigen::Matrix3f R = T.block<3, 3>(0, 0);
    const double t_err = T.block<3, 1>(0, 3).cast<double>().norm();
    const double r_err = std::fabs(static_cast<double>(Eigen::AngleAxisf(R).angle()));
    std::printf("用例1: converged=%d iters=%d t_err=%.6f r_err=%.6f\n", ndt.hasConverged() ? 1 : 0,
                ndt.getFinalNumIteration(), t_err, r_err);
    Check(ndt.hasConverged(), "用例1 端到端：hasConverged 应为真");
    Check(AllFinite(T), "用例1 端到端：最终变换应有限");
    Check(t_err < 0.01, "用例1 端到端：平移误差应小于 0.01 m（目标与源为同一点集）");
    Check(r_err < 0.005, "用例1 端到端：旋转误差应小于 0.005 rad");
  }

  // 用例 1b：单位变换为驻点的依据——对称点簇下梯度应为零（相对 Hessian 规模）
  {
    Vec6 x_id = Vec6::Zero();
    Vec6 g;
    Mat6 H;
    const double s = Derivs(ndt, *source, x_id, g, H, true);
    std::printf("用例1b: score=%.6f |g|=%.6e |H|=%.6e\n", s, g.norm(), H.norm());
    Check(std::isfinite(s) && H.norm() > 1e-9, "用例1b 驻点：单位变换处得分有限、Hessian 非零");
    Check(g.norm() <= 1e-5 * H.norm(), "用例1b 驻点：单位变换处梯度应接近零（对称点簇，梯度范数 ≤ 1e-5 × Hessian 范数，仅为相对尺度归一化，不等价于牛顿步上界）");
  }

  // 用例 2/3：历史策略保留过大的初始试探步；得分、Hessian 和输出点云必须对应实际选用的位置。
  {
    Vec6 x0;
    x0 << 0.15, 0.0, 0.0, 0.0, 0.0, 0.0;  // 源沿 x 偏移 0.15，真实最优在 x=0
    Vec6 g0;
    Mat6 H0;
    const double score0 = Derivs(ndt, *source, x0, g0, H0, true);

    Vec6 dir;
    dir << -1.0, 0.0, 0.0, 0.0, 0.0, 0.0;  // 指向最优的方向；若函数判为非下降方向会自行反向（dir 为引用）
    const double step_init = 0.6, step_max = 0.6, step_min = 1e-3;  // 0.6 远大于最优步长 0.15 → 明显过冲

    double score = score0;
    Vec6 g = g0;
    Mat6 H = H0;
    CloudT trans;
    const double a_t = ndt.computeStepLengthMT(x0, dir, step_init, step_max, step_min, score, g, H, trans);
    std::printf("用例2: score0=%.6f score_t=%.6f a_t=%.6f |H|=%.6e\n", score0, score, a_t, H.norm());

    Check(std::isfinite(a_t) && a_t > 0.0, "用例2 线搜索：返回步长应为有限正数");
    Check(a_t == step_init,
          "用例2 历史步长：正常区间保留截断后的初始牛顿步，不执行内层细化");
    Check(score < score0, "用例2 场景：所选步长确实过冲，未暗中恢复改善得分的线搜索");

    Vec6 x_fin = x0 + dir * a_t;
    Vec6 g2;
    Mat6 H2;
    const double score2 = Derivs(ndt, *source, x_fin, g2, H2, true);
    const double h_scale = std::max(H2.norm(), 1e-12);
    Check(H.norm() > 1e-9,
          "用例3 Hessian：所选初始试探步的 Hessian 不应为零");
    Check((H - H2).norm() <= 1e-8 * (1.0 + h_scale), "用例3 Hessian：返回的 Hessian 应与最终位置重算结果一致");
    Check(std::fabs(score - score2) <= 1e-9 * (1.0 + std::fabs(score2)), "用例3 一致性：返回得分应与最终位置重算一致");
    Check((g - g2).norm() <= 1e-9 * (1.0 + g2.norm()), "用例3 一致性：返回梯度应与最终位置重算一致");
    Check(trans.size() == source->size(), "用例3 一致性：trans_cloud 点数应与输入源一致");
  }

  // 用例 4：梯度中心差分核对解析 Hessian。先用实际 round 键证明所有评估位姿下无邻域跳变，再做差分。
  {
    Vec6 xc;
    xc << 0.006, -0.005, 0.004, 0.002, -0.0025, 0.003;  // 角度绝对值差分后仍 ≥1e-3，远离小角近似阈值 1e-4
    const std::vector<Eigen::Vector3i> keys_id = KeysOf(*source);
    const std::vector<Eigen::Vector3i> keys_c = KeysOf(Transformed(*source, xc));
    bool keys_ok = Check(SameKeys(keys_id, keys_c), "用例4 前提：基准位姿下所有点的 round 体素键与单位变换一致");

    Vec6 g0;
    Mat6 H0;
    const double s0 = Derivs(ndt, *source, xc, g0, H0, true);
    Mat6 Hfd = Mat6::Zero();
    Mat6 no_hess = Mat6::Ones();
    for (int i = 0; i < 6 && keys_ok; ++i) {
      const double h = (i < 3) ? 2e-3 : 1e-3;
      Vec6 xp = xc, xm = xc;
      xp(i) += h;
      xm(i) -= h;
      keys_ok = Check(SameKeys(keys_c, KeysOf(Transformed(*source, xp))) &&
                          SameKeys(keys_c, KeysOf(Transformed(*source, xm))),
                      "用例4 前提：差分位姿 " + std::to_string(i) + " 下所有点的 round 体素键不变（无邻域跳变）");
      if (!keys_ok) break;
      Vec6 gp, gm;
      Mat6 Hp, Hm;
      Derivs(ndt, *source, xp, gp, Hp, false);
      Derivs(ndt, *source, xm, gm, Hm, false);
      Hfd.col(i) = (gp - gm) / (2.0 * h);
      if (i == 0) no_hess = Hp;
    }
    if (keys_ok) {
      const double rel = (Hfd - H0).norm() / std::max(H0.norm(), 1e-12);
      std::printf("用例4: score=%.6f |H|=%.6e |Hfd-H|/|H|=%.6f\n", s0, H0.norm(), rel);
      Check(std::isfinite(s0) && H0.norm() > 1e-9, "用例4 Hessian：解析 Hessian 应有限且非零");
      Check(rel < 0.05,
            "用例4 Hessian：无邻域跳变、远离边界的平滑区间内，梯度中心差分与解析 Hessian 相对偏差应小于 5%");
      Check(no_hess.norm() == 0.0, "用例4 前提：compute_hessian=false 时返回零 Hessian（修复所依据的事实）");
    } else {
      std::printf("用例4: 因 round 键前提不成立而中止差分核对，场景需重新设计\n");
    }
  }

  // 用例 5：迭代预算耗尽不等于收敛（修复不得改变该语义）
  {
    ndt.setMaximumIterations(1);
    ndt.setTransformationEpsilon(1e-12);
    Vec6 guess;
    guess << 0.4, 0.3, 0.0, 0.0, 0.0, 0.1;
    const Eigen::Matrix4f T = RunAlign(ndt, guess);
    std::printf("用例5: converged=%d iters=%d\n", ndt.hasConverged() ? 1 : 0, ndt.getFinalNumIteration());
    Check(ndt.getFinalNumIteration() == 1, "用例5 预算：应恰好用尽迭代预算（上限 1，实际迭代 1 次），否则未收敛可能来自其他提前退出路径，本用例无效");
    Check(!ndt.hasConverged(), "用例5 预算：达到迭代预算且步长未低于容差时 hasConverged 不应为真");
    Check(AllFinite(T), "用例5 预算：预算耗尽时最终变换仍应有限");
    ndt.setMaximumIterations(50);
    ndt.setTransformationEpsilon(1e-4);
  }

  if (g_failures != 0) {
    std::printf("结果：共 %d 项检查失败\n", g_failures);
    return 1;
  }
  std::printf("结果：全部检查通过（仅代表合成用例通过，不代表原始事件已解决）\n");
  return 0;
}
