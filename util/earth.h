#ifndef EARTH_H
#define EARTH_H

/* WGS84椭球模型参数
   NOTE:如果使用其他椭球模型需要修改椭球参数 */
const double WGS84_WIE = 7.2921151467E-5; /* 地球自转角速度*/
const double WGS84_F = 0.0033528106647474805; /* 扁率 */
const double WGS84_RA = 6378137.0000000000; /* 长半轴a */
const double WGS84_RB = 6356752.3142451793; /* 短半轴b */
const double WGS84_GM0 = 398600441800000.00; /* 地球引力常数 */
const double WGS84_E1 = 0.0066943799901413156; /* 第一偏心率平方 */
const double WGS84_E2 = 0.0067394967422764341; /* 第二偏心率平方 */

namespace Earth
{
    inline Eigen::Vector2d meridianPrimeVerticalRadius(double lat)
    {
        double tmp, sqrttmp;

        tmp = sin(lat);
        tmp *= tmp;
        tmp = 1 - WGS84_E1 * tmp;
        sqrttmp = sqrt(tmp);

        return {WGS84_RA * (1 - WGS84_E1) / (sqrttmp * tmp), WGS84_RA / sqrttmp};
    }

    inline double gravity(const Eigen::Vector3d &blh)
    {
        double sinphi = sin(blh[0]);
        double sin2 = sinphi * sinphi;
        double sin4 = sin2 * sin2;

        // normal gravity at equator, 赤道处正常重力
        double gamma_a = 9.7803267715;
        // series expansion of normal gravity at given latitude, 给定纬度处正常重力的级数展开
        double gamma_0 = gamma_a * (1 + 0.0052790414 * sin2 + 0.0000232718 * sin4 + 0.0000001262 * sin2 * sin4 +
                                    0.0000000007 * sin4 * sin4);
        // changes of normal gravity with height, 正常重力随高度变化
        double gamma = gamma_0 - (3.0877e-6 - 4.3e-9 * sin2) * blh[2] + 0.72e-12 * blh[2] * blh[2];

        return gamma;
    }

    /* n系(导航坐标系)到e系(地心地固坐标系)转换四元数 */
    inline Eigen::Quaterniond qne(const Eigen::Vector3d &blh)
    {
        Eigen::Quaterniond quat;

        double coslon, sinlon, coslat, sinlat;

        coslon = cos(blh[1] * 0.5);
        sinlon = sin(blh[1] * 0.5);
        coslat = cos(-M_PI * 0.25 - blh[0] * 0.5);
        sinlat = sin(-M_PI * 0.25 - blh[0] * 0.5);

        quat.w() = coslat * coslon;
        quat.x() = -sinlat * sinlon;
        quat.y() = sinlat * coslon;
        quat.z() = coslat * sinlon;

        return quat.normalized();
    }

    /* 从n系到e系转换四元数得到纬度和经度 */
    inline  Eigen::Vector3d blh(const Eigen::Quaterniond &qne, double height) {
        return {-2 * atan(qne.y() / qne.w()) - M_PI * 0.5, 2 * atan2(qne.z(), qne.w()), height};
    }

    /* n系相对位置转大地坐标相对位置 */
    inline Eigen::Matrix3d DRi(const Eigen::Vector3d &blh) {
        Eigen::Matrix3d dri = Eigen::Matrix3d::Zero();

        Eigen::Vector2d rmn = meridianPrimeVerticalRadius(blh[0]);

        dri(0, 0) = 1.0 / (rmn[0] + blh[2]);
        dri(1, 1) = 1.0 / ((rmn[1] + blh[2]) * cos(blh[0]));
        dri(2, 2) = -1;
        return dri;
    }

};

#endif
