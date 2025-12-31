#ifndef MEMS_INS_GNSS_ESKF_ROTATION_H
#define MEMS_INS_GNSS_ESKF_ROTATION_H

#include <iostream>
#include <eigen/Eigen/Dense>

namespace Rotation
{
    inline Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d &vector)
    {
        Eigen::Matrix3d mat;
        mat << 0, -vector(2), vector(1), vector(2), 0, -vector(0), -vector(1), vector(0), 0;
        return mat;
    }

    inline Eigen::Quaterniond rotvec2quaternion(const Eigen::Vector3d &rotvec)
    {
        double angle = rotvec.norm();
        Eigen::Vector3d vec = rotvec.normalized();
        return Eigen::Quaterniond(Eigen::AngleAxisd(angle, vec));
    }

    inline Eigen::Vector3d quaternion2vector(const Eigen::Quaterniond &quaternion)
    {
        Eigen::AngleAxisd axisd(quaternion);
        return axisd.angle() * axisd.axis();
    }

    inline Eigen::Matrix3d quaternion2matrix(const Eigen::Quaterniond &quaternion)
    {
        return quaternion.toRotationMatrix();
    }

    inline Eigen::Vector3d matrix2euler(const Eigen::Matrix3d &dcm)
    {
        Eigen::Vector3d euler;

        euler[1] = atan(-dcm(2, 0) / sqrt(dcm(2, 1) * dcm(2, 1) + dcm(2, 2) * dcm(2, 2)));

        if (dcm(2, 0) <= -0.999)
        {
            euler[0] = 0;
            euler[2] = atan2((dcm(1, 2) - dcm(0, 1)), (dcm(0, 2) + dcm(1, 1)));
            std::cout << "[WARNING] Rotation::matrix2euler: Singular Euler Angle! Set the roll angle to 0!" <<
                    std::endl;
        }
        else if (dcm(2, 0) >= 0.999)
        {
            euler[0] = 0;
            euler[2] = M_PI + atan2((dcm(1, 2) + dcm(0, 1)), (dcm(0, 2) - dcm(1, 1)));
            std::cout << "[WARNING] Rotation::matrix2euler: Singular Euler Angle! Set the roll angle to 0!" <<
                    std::endl;
        }
        else
        {
            euler[0] = atan2(dcm(2, 1), dcm(2, 2));
            euler[2] = atan2(dcm(1, 0), dcm(0, 0));
        }

        // heading 0~2PI
        if (euler[2] < 0)
        {
            euler[2] = M_PI * 2 + euler[2];
        }

        return euler;
    }
}
#endif //MEMS_INS_GNSS_ESKF_ROTATION_H
