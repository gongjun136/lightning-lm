//
// Created by xiang on 25-4-14.
//

#ifndef FASTER_LIO_POSE6D_H
#define FASTER_LIO_POSE6D_H

#include "common/eigen_types.h"

namespace lightning {

struct Pose6D {
    Pose6D() = default;

    Pose6D(const double t, const Vec3d &a, const Vec3d &g, const Vec3d &v, const Vec3d &p, const Mat3d &R) {
        offset_time = t;
        acc = a;
        gyr = g;
        vel = v;
        pos = p;
        rot = R;
    };

    double offset_time = 0;         // 相对于点云开始时间的时间偏移量，用于点云运动畸变校正
    Vec3d acc = Vec3d::Zero();      // 加速度（世界坐标系下，区间平均值，经过零偏校正并转换到世界坐标系，并减去重力）
    Vec3d gyr = Vec3d::Zero();      // 角速度（IMU坐标系下，区间平均值，经过零偏校正）
    Vec3d vel = Vec3d::Zero();      // 速度（世界坐标系下，区间起始时刻的状态）
    Vec3d pos = Vec3d::Zero();      // 位置（世界坐标系下，区间起始时刻的状态）
    Mat3d rot = Mat3d::Identity();  // 旋转矩阵（从IMU坐标系到世界坐标系，区间起始时刻的状态）
};

}  // namespace lightning

#endif  // FASTER_LIO_POSE6D_H
