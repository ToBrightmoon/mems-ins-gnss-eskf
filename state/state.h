#ifndef STATE_H
#define STATE_H

#include <eigen/Eigen/Dense>

namespace State
{
    struct Attitude
    {
        Eigen::Quaterniond qbn;
        Eigen::Matrix3d cbn;
        Eigen::Vector3d euler;
    };

    struct ImuError
    {
        Eigen::Vector3d gyrBias;
        Eigen::Vector3d accBias;
    };

    struct Pva
    {
        Eigen::Vector3d pos;
        Eigen::Vector3d vel;
        Attitude att;
    };

    struct  NormalState
    {
        Pva pva;
        ImuError imuErr;
    };

}
#endif //STATE_H
