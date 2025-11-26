//
// Created by xiang on 2022/2/15.
//
#pragma once

#include "core/lightning_math.hpp"

namespace lightning {

/**
 * sphere 2
 * 本身是三维的，更新量是二维的一种东东
 * 底下那堆计算我反正看不懂，别问我是啥，问就是照抄的
 *
 * @todo S2的长度应当可以外部指定
 */
struct S2 {
    using SO3 = Sophus::SO3d;

    static constexpr int den_ = 98090;
    static constexpr int num_ = 10000;
    static constexpr double length_ = double(den_) / double(num_);
    Vec3d vec_;

   public:
    S2() { vec_ = length_ * Vec3d(1, 0, 0); }
    S2(const Vec3d &vec) : vec_(vec) {
        vec_.normalize();
        vec_ = vec_ * length_;
    }

    double operator[](int idx) const { return vec_[idx]; }

    Eigen::Matrix<double, 3, 3> S2_hat() const {
        Eigen::Matrix<double, 3, 3> skew_vec;
        skew_vec << double(0), -vec_[2], vec_[1], vec_[2], double(0), -vec_[0], -vec_[1], vec_[0], double(0);
        return skew_vec;
    }

    /**
     * S2_Mx: 计算S2流形上boxplus操作的雅可比矩阵 ∂(x ⊕ δ)/∂δ
     * 用于ESKF中重力向量误差的线性化传播
     */
    Eigen::Matrix<double, 3, 2> S2_Mx(const Eigen::Matrix<double, 2, 1> &delta) const {
        Eigen::Matrix<double, 3, 2> res;
        Eigen::Matrix<double, 3, 2> Bx = S2_Bx();  // 获取S2切空间的基矩阵

        // 小角度近似：使用一阶泰勒展开
        if (delta.norm() < 1e-5) {
            res = -SO3::hat(vec_) * Bx;  // 线性近似下的雅可比
        } else {
            // 大角度：使用完整的李群雅可比计算
            Vec3d Bu = Bx * delta;                // 将2D增量映射到3D切空间
            SO3 exp_delta = math::exp(Bu, 0.5f);  // 指数映射到SO3
            /**
             * Derivation of d(Exp(Bx dx)x)/d(dx)=d(Exp(Bu)x)/d(dx):
             *  d(Exp(Bu)x)/d(dx)=d(Exp(Bu)x)/d(Bu) Bx; then
             *  d(Exp(Bu)x)/d(Bu)=d[Exp(Bu+dBu)x]/d(dBu)=d[Exp(Jl(Bu)dBu)Exp(Bu)x]/d(dBu)=d[(Jl(Bu)dBu)^Exp(Bu)x]/d(dBu)
             *   =d[-(Exp(Bu)x)^Jl(Bu)dBu]/d(dBu)=-Exp(Bu)x^Exp(-Bu)Jl(Bu);
             *    for Exp(x+dx)=Exp(x)Exp(Jr(x)dx)=Exp(Jl(x)dx)Exp(x)=Exp(x)Exp(Exp(-x)Jl(x)dx) =>
             *    Exp(-x)Jl(x)=Jr(x)=Jl(-x) =>
             *   =-Exp(Bu)x^Jl(-Bu) => d(Exp(Bu)x)/d(dx)= -Exp(Bu) x^ Jl(Bu)^T Bx or A_matrix is just Jl()
             */
            res = -exp_delta.matrix() * SO3::hat(vec_) * math::A_matrix(Bu).transpose() * Bx;
        }
        return res;
    }

    /**
     * Bx两个列向量为正切空间的局部坐标系
     *
     * S2 = [a b c], l = length
     * Bx = [
     * -b              -c
     * l-bb/(l+a)      -bc/(l+a)
     * -bc/(l+a)        l-cc/(l+a)
     * ] / l
     * Derivation of origin MTK: (112) Rv = [v Bx] = [x -r 0; y xc -s; z xs c]
     *  where c=cos(alpha),s=sin(alpha),r=sqrt(y^2+z^2),||v||=1:
     *  For y-axis or (0,1,0) in local frame is mapped to (-r;xc;xs) in world frame,
     *  meaning x/OG(or gravity vector) projects to A of yz plane,
     *  and this projecting line will intersect the perpendicular plane of OG at O at the point B,
     *  then from OB we can found OC with ||OC||=1, which is just the y-axis, then we get its coordinate in world frame:
     *  ∠AOG+∠AGO=pi/2=∠AOG+∠AOB => ∠AOB=∠AGO, for sin(∠AGO)=r/||v||=r => ||OC||*sin(∠AOB)=r;||OC||*cos(∠AOB)=x =>
     *  (-r;xc;xs) is the y-axis coordinate in world frame, then z-axis is easy to get
     * Derivation of current MTK with S2 = [x y z], ||S2|| = 1:
     *  just a rotation of origin one around x/OG axis to meet y-axis coordinate in word frame to be (-y;a;c),
     *  then z-axis must be (+-z;b;d) for x^+y^+z^=1, where a,b,c,d is just variable to be solved
     *  for current MTK chooses clockwise rotation(meaning -z):
     *  Rv=[x -y -z;
     *      y  a  b;
     *      z  c  d]
     *  then for -xy+ya+zc=0=xy-ya-zb => b=c
     *  for yz+ab+cd=0; b=c => a=-d-yz/c; for -xz+yb+zd=0; b=c => d=x-yc/z
     *  for -xy+ya+zc=0; a=-d-yz/c=yc/z-x-yz/c => -xy+y^2c/z-xy-y^2z/c+zc=0 => (z+y^2/z)c^2 -2xy c -y^2z = 0 =>
     *  c=[2xy +- sqrt(4x^2y^2 + 4(z+y^2/z)y^2z)]/[2(z+y^2/z)]; for z^2+y^2=1-x^2 =>
     *  c=[xy +- y]z/(1-x^2), for rotation is clocewise, thus c<0 => c=(xy-y)z/(1-x^2)=-yz/(1+x)
     *  then b,a,d is easy to get and also if ||S2||=l, it is easy to prove c=-ylzl/(l+lx)/l=-bc/(l+a)/l
     *
     * Bx: 计算S2流形在当前点的切空间基矩阵 (3×2)
     *
     * 数学原理：
     * - S2流形上的每个点都有一个2维切空间（单位球面的切平面）
     * - Bx的两个列向量构成这个切空间的正交基
     * - 用于将2维误差增量映射到3维空间中进行李群操作
     *
     * 几何解释：
     * - 对于重力向量 g = [x, y, z]^T，长度为 l
     * - 切平面与重力向量垂直，通过坐标原点
     * - Bx将2维局部坐标系增量映射到3D世界坐标系
     */
    Eigen::Matrix<double, 3, 2> S2_Bx() const {
        Eigen::Matrix<double, 3, 2> res;
        if (vec_[0] + length_ > 1e-5) {
            res << -vec_[1], -vec_[2], length_ - vec_[1] * vec_[1] / (length_ + vec_[0]),
                -vec_[2] * vec_[1] / (length_ + vec_[0]), -vec_[2] * vec_[1] / (length_ + vec_[0]),
                length_ - vec_[2] * vec_[2] / (length_ + vec_[0]);
            res /= length_;
        } else {
            res = Eigen::Matrix<double, 3, 2>::Zero();
            res(1, 1) = -1;
            res(2, 0) = 1;
        }
        return res;
    }

    /**
     * S2_Nx_yy: 计算S2流形boxminus操作的雅可比矩阵 ∂(x ⊖ y)/∂x，当x=y时的简化版本
     * 用于ESKF中误差状态的线性化和协方差传播
     *
     * 数学含义：
     * - 当两个S2向量相等时，计算boxminus对第一个向量的偏导数
     * - 结果为2×3矩阵，将3D空间误差映射到2D切空间
     * - 长度因子1/length²确保单位球面约束
     */
    Eigen::Matrix<double, 2, 3> S2_Nx_yy() const {
        Eigen::Matrix<double, 2, 3> res;
        Eigen::Matrix<double, 3, 2> Bx = S2_Bx();                       // 获取切空间基矩阵
        res = 1 / length_ / length_ * Bx.transpose() * SO3::hat(vec_);  // 雅可比 = (1/|g|²) * Bx^T * [g]×
        return res;
    }

    /**
     * oplus: S2流形上的广义加法运算 (李群右乘)
     * vec_ = vec_ ⊕ delta = Exp(delta*0.5) * vec_
     *
     * 数学含义：
     * - 将3D增量向量delta通过SO3指数映射施加到当前向量
     * - 使用scale参数控制增量的大小
     * - 结果自动保持S2流形约束（固定长度）
     *
     * 在ESKF中用于：名义状态更新 = 名义状态 ⊕ 误差状态
     */
    void oplus(const Vec3d &delta, double scale = 1.0) { vec_ = math::exp(delta, scale * 0.5) * vec_; }

    /**
     * boxminus: S2流形上的广义减法运算
     * delta = vec_ ⊖ other = Log(other^T * vec_)
     * 返回两个S2向量之间的最小旋转弧长在切空间的表示
     *
     * 数学原理：
     * - 计算从other到vec_的最小旋转
     * 使用atan2(v_sin, v_cos)得到旋转角度theta
     * - 通过Bx矩阵将旋转映射到2维切空间
     * - 结果为2D向量，表示在other切空间的误差
     *
     * 在ESKF中用于：误差状态 = 名义状态 ⊖ 先验状态
     */
    Vec2d boxminus(const S2 &other) const {
        Vec2d res;
        // 计算两个向量的夹角：sin(θ) = |[vec_]×other|, cos(θ) = vec_^T * other
        double v_sin = (SO3::hat(vec_) * other.vec_).norm();
        double v_cos = vec_.transpose() * other.vec_;
        double theta = std::atan2(v_sin, v_cos);  // 反正切得到角度

        // 处理数值稳定性：小角度和大角度的特殊情况
        if (v_sin < 1e-5) {
            if (std::fabs(theta) > 1e-5) {
                res[0] = 3.1415926;  // π: 完全相反
                res[1] = 0;
            } else {
                res[0] = 0;  // 0: 完全相同
                res[1] = 0;
            }
        } else {
            // 大角度：使用完整公式计算切空间误差
            S2 other_copy = other;
            Eigen::Matrix<double, 3, 2> Bx = other_copy.S2_Bx();
            res = theta / v_sin * Bx.transpose() * SO3::hat(other.vec_) * vec_;
        }
        return res;
    }

    /**
     * boxplus: S2流形上的广义加法运算（2维版本）
     * vec_ = vec_ ⊕ delta = Exp(Bx * delta * scale/2) * vec_
     *
     * 与oplus的区别：
     * - oplus接受3D增量，自动处理切空间映射
     * boxplus接受2维增量，需要显式使用Bx矩阵
     *
     * 数学含义：
     * - Bx将2维增量映射到3D切空间向量
     * - SO3指数映射将3D向量转换为旋转矩阵
     * - 旋转矩阵作用在当前向量上实现更新
     * - 结果自动保持S2流形约束
     *
     * 在ESKF中用于：先验状态更新 = 名义状态 ⊕ 误差增量
     */
    void boxplus(const Vec2d &delta, double scale = 1) {
        Eigen::Matrix<double, 3, 2> Bx = S2_Bx();    // 获取切空间基矩阵
        SO3 res = math::exp(Bx * delta, scale / 2);  // 将2D增量映射为SO3旋转
        vec_ = res.matrix() * vec_;                  // 应用旋转更新向量
    }
};

}  // namespace lightning
