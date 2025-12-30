#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <eigen/Eigen/Geometry>

namespace Measurement
{
    using V3d = Eigen::Vector3d;
    struct IMU
    {
        double time;
        double dt;
        V3d thetaD;
        V3d velD;
    };
}
#endif //MEASUREMENT_H
