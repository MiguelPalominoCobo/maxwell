#include "Hopfion.h"
#include "Types.h"
#include <iostream>

Hopfion::Hopfion(std::size_t pIn, std::size_t qIn)
{
    p = pIn;
    q = qIn;
}

Hopfion::FieldEH Hopfion::evaluate(double time, Vec3 position) const
{
    double x = position[0];
    double y = position[1];
    double z = position[2];

    Hopfion::Vec3 E;
    Hopfion::Vec3 H;

    H[0] = fieldHx(time, x, y, z);
    H[1] = fieldHy(time, x, y, z);
    H[2] = fieldHz(time, x, y, z);
    E[0] = fieldEx(time, x, y, z);
    E[1] = fieldEy(time, x, y, z);
    E[2] = fieldEz(time, x, y, z);

    // Hopfion::FieldEH CamposEH = <E, H>;
    // Hopfion::FieldEH CamposEH(E, H);
    Hopfion::FieldEH CamposEH = std::make_pair(E, H);

    return CamposEH;

}
//
//double Hopfion::evalEX(const double time, const double x, const double y, const double z) const
//{
//    return fieldEx(time, x, y, z);
//}
//
//double Hopfion::evalEY(const double time, const double x, const double y, const double z) const
//{
//    return Hopfion::fieldEy(time, x, y, z);
//}
//
//double Hopfion::evalEZ(const double time, const double x, const double y, const double z) const
//{
//    return Hopfion::fieldEz(time, x, y, z);
//}
//
//double Hopfion::evalHX(const double time, const double x, const double y, const double z) const
//{
//    return Hopfion::fieldHx(time, x, y, z);
//}
//
//double Hopfion::evalHY(const double time, const double x, const double y, const double z) const
//{
//    return Hopfion::fieldHy(time, x, y, z);
//}
//
//double Hopfion::evalHZ(const double time, const double x, const double y, const double z) const
//{
//    return Hopfion::fieldHz(time, x, y, z);
//}