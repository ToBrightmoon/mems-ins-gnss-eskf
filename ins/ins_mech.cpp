#include "ins_mech.h"

#include "util/earth.h"
#include "util/rotation.h"

namespace Ins
{
    static Eigen::Vector3d posUpdate(const State::Pva &prePva,const Eigen::Vector3d &currVel,
                                     const Measurement::IMU &currImu)
    {
        Eigen::Vector3d temp1, temp2, midvel, midpos;
        Eigen::Quaterniond qne, qee, qnn;

        // 重新计算中间时刻的速度和位置
        midvel = (currVel + prePva.vel) / 2;
        midpos = prePva.pos + Earth::DRi(prePva.pos) * midvel * currImu.dt / 2;

        // 重新计算中间时刻地理参数
        Eigen::Vector2d rmrn;
        Eigen::Vector3d wie_n, wen_n;
        rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
        wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
        wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
            -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

        // 重新计算 k时刻到k-1时刻 n系旋转矢量
        temp1 = (wie_n + wen_n) * currImu.dt;
        qnn   = Rotation::rotvec2quaternion(temp1);
        // e系转动等效旋转矢量 (k-1时刻k时刻，所以取负号)
        temp2 << 0, 0, -WGS84_WIE * currImu.dt;
        qee = Rotation::rotvec2quaternion(temp2);

        // 位置更新完成
        qne           = Earth::qne(prePva.pos);
        qne           = (qee * qne * qnn).normalized();
        auto currHeigh = prePva.pos[2] - midvel[2] * currImu.dt;
        return Earth::blh(qne, currHeigh);
    };

    static Eigen::Vector3d velUpdate(const State::Pva &prePva, const Measurement::IMU &preImu,
                                     const Measurement::IMU &currImu)
    {
        Eigen::Vector3d d_vfb, d_vfn, d_vgn, gl, midvel, midpos;
        Eigen::Vector3d temp1, temp2, temp3;
        Eigen::Matrix3d cnn, I33 = Eigen::Matrix3d::Identity();
        Eigen::Quaterniond qne, qee, qnn, qbb, q1, q2;

        // 计算地理参数，子午圈半径和卯酉圈半径，地球自转角速度投影到n系, n系相对于e系转动角速度投影到n系，重力值
        Eigen::Vector2d rmrn = Earth::meridianPrimeVerticalRadius(prePva.pos(0));
        Eigen::Vector3d wie_n, wen_n;
        wie_n << WGS84_WIE * cos(prePva.pos[0]), 0, -WGS84_WIE * sin(prePva.pos[0]);
        wen_n << prePva.vel[1] / (rmrn[1] + prePva.pos[2]), -prePva.vel[0] / (rmrn[0] + prePva.pos[2]),
                -prePva.vel[1] * tan(prePva.pos[0]) / (rmrn[1] + prePva.pos[2]);
        double gravity = Earth::gravity(prePva.pos);

        // 旋转效应和双子样划桨效应
        temp1 = currImu.thetaD.cross(currImu.velD) / 2;
        temp2 = preImu.thetaD.cross(currImu.velD) / 12;
        temp3 = preImu.velD.cross(currImu.thetaD) / 12;

        // b系比力积分项
        d_vfb = currImu.velD + temp1 + temp2 + temp3;

        // 比力积分项投影到n系
        temp1 = (wie_n + wen_n) * currImu.dt / 2;
        cnn = I33 - Rotation::skewSymmetric(temp1);
        d_vfn = cnn * prePva.att.cbn * d_vfb;

        // 计算重力/哥式积分项
        gl << 0, 0, gravity;
        d_vgn = (gl - (2 * wie_n + wen_n).cross(prePva.vel)) * currImu.dt;

        // 得到中间时刻速度
        // velocity at k-1/2
        midvel = prePva.vel + (d_vfn + d_vgn) / 2;

        // 外推得到中间时刻位置
        qnn = Rotation::rotvec2quaternion(temp1);
        temp2 << 0, 0, -WGS84_WIE * currImu.dt / 2;
        qee = Rotation::rotvec2quaternion(temp2);
        qne = Earth::qne(prePva.pos);
        qne = (qee * qne * qnn).normalized();
        midpos[2] = prePva.pos[2] - midvel[2] * currImu.dt / 2;
        midpos = Earth::blh(qne, midpos[2]);

        // 重新计算中间时刻的rmrn, wie_e, wen_n
        rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
        wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
        wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
                -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

        // 重新计算n系下平均比力积分项
        temp3 = (wie_n + wen_n) * currImu.dt / 2;
        cnn = I33 - Rotation::skewSymmetric(temp3);
        d_vfn = cnn * prePva.att.cbn * d_vfb;

        // 重新计算重力、哥式积分项
        gl << 0, 0, Earth::gravity(midpos);
        d_vgn = (gl - (2 * wie_n + wen_n).cross(midvel)) * currImu.dt;

        return prePva.vel + d_vfn + d_vgn;
    };

    static Eigen::Quaterniond attUpdate(const State::Pva &prePva, const Measurement::IMU &preImu,
                                        const Measurement::IMU &currImu,const Eigen::Vector3d &currPos, const Eigen::Vector3d &currVel)
    {
        Eigen::Quaterniond qne_pre, qne_cur, qne_mid, qnn, qbb;
        Eigen::Vector3d temp1, midpos, midvel;

        // 重新计算中间时刻的速度和位置
        midvel = (prePva.vel + currVel) / 2;
        qne_pre   = Earth::qne(prePva.pos);
        qne_cur   = Earth::qne(currPos);
        temp1     = Rotation::quaternion2vector((qne_cur.inverse() * qne_pre).normalized());
        qne_mid   = (qne_pre * Rotation::rotvec2quaternion(temp1 / 2).inverse()).normalized();
        midpos[2] = (currPos[2] + prePva.pos[2]) / 2;
        midpos    = Earth::blh(qne_mid, midpos[2]);

        // 重新计算中间时刻地理参数
        Eigen::Vector2d rmrn;
        Eigen::Vector3d wie_n, wen_n;
        rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
        wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
        wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
            -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

        // 计算n系的旋转四元数 k-1时刻到k时刻变换
        temp1 = -(wie_n + wen_n) * currImu.dt;
        qnn   = Rotation::rotvec2quaternion(temp1);

        // 计算b系旋转四元数 补偿二阶圆锥误差
        temp1 = currImu.thetaD + preImu.thetaD.cross(currImu.thetaD) / 12;
        qbb   = Rotation::rotvec2quaternion(temp1);

        // 姿态更新完成
        return  (qnn * prePva.att.qbn * qbb).normalized();
    };

    State::Pva insMech(const State::Pva &prePva, const Measurement::IMU &preImu, const Measurement::IMU &currImu)
    {
        State::Pva currPva;
        currPva.vel = velUpdate(prePva,preImu,currImu);
        currPva.pos = posUpdate(prePva,currPva.vel,currImu);
        currPva.att.qbn = attUpdate(prePva,preImu,currImu,currPva.pos,currPva.vel);
        currPva.att.cbn = Rotation::quaternion2matrix(currPva.att.qbn);
        currPva.att.euler = Rotation::matrix2euler(currPva.att.cbn);
        return currPva;
    };
};
