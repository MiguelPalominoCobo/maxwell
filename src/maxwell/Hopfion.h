#pragma once

#include <functional>
#include <array>
#include <complex>
#include "Types.h"

#include "mfem.hpp"


using namespace maxwell;

class Hopfion {
public:

	typedef std::array<double, 3> Vec3;
	typedef std::pair<Vec3, Vec3> FieldEH;

	Hopfion(std::size_t p, std::size_t q);

	FieldEH evaluate(double time, Vec3 pos) const;

    double getEX(const double time, const double x, const double y, const double z) const { return fieldEx(time, x, y, z); }
    double getEY(const double time, const double x, const double y, const double z) const { return fieldEy(time, x, y, z); }
    double getEZ(const double time, const double x, const double y, const double z) const { return fieldEz(time, x, y, z); }
    double getHX(const double time, const double x, const double y, const double z) const { return fieldHx(time, x, y, z); }
    double getHY(const double time, const double x, const double y, const double z) const { return fieldHy(time, x, y, z); }
    double getHZ(const double time, const double x, const double y, const double z) const { return fieldHz(time, x, y, z); }



private:

	std::size_t p, q;


    double fieldEx(const double t, const double x, const double y, const double z) const
    {
        double Ex;
        std::complex<double> I(0, 1);
        std::complex<double> Fx;

        Fx = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) * (-4 * x * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2) + 2.0 * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * x / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2)) - 2 * x * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2));
        Ex = real(Fx);

        return Ex;
    }

    double fieldHx(const double t, const double x, const double y, const double z) const
    {
        double Hx;
        std::complex<double> I(0, 1);
        std::complex<double> Fx;

        Fx = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) * (-4 * x * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2) + 2.0 * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * x / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2)) - 2 * x * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow((-I + t), 2), 2));
        Hx = imag(Fx);

        return Hx;
    }


    double fieldEy(const double t, const double x, const double y, const double z) const
    {
        double Ey;
        std::complex<double> I(0, 1);
        std::complex<double> Fy;

        Fy = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) * (-4 * y * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 2.0 * I * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * y / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2 * y * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Ey = real(Fy);

        return Ey;
    }

    double fieldHy(const double t, const double x, const double y, const double z) const
    {
        double Hy;
        std::complex<double> I(0, 1);
        std::complex<double> Fy;

        Fy = I * ((-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) * (-4 * y * (x - I * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 2.0 * I * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), -1)) - 4.0 * (x - I * y) * (-I + t) * (2 * y / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2 * y * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Hy = imag(Fy);

        return Hy;
    }


    double fieldEz(const double t, const double x, const double y, const double z) const
    {
        double Ez;
        std::complex<double> I(0, 1);
        std::complex<double> Fz;

        Fz = I * (-4 * z * (x - I * y) * (-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 4.0 * (x - I * y) * (-I + t) * ((2.0 * I + 2.0 * z) / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2.0 * z * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Ez = real(Fz);

        return Ez;
    }

    double fieldHz(const double t, const double x, const double y, const double z) const
    {
        double Hz;
        std::complex<double> I(0, 1);
        std::complex<double> Fz;

        Fz = I * (-4 * z * (x - I * y) * (-2 * t / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) + 2.0 * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) * (-I + t) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2) - 4.0 * (x - I * y) * (-I + t) * ((2.0 * I + 2.0 * z) / (pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2)) - 2.0 * z * (-1.0 + 2.0 * I * z - pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - pow(-I + t, 2), 2));
        Hz = imag(Fz);

        return Hz;
    }

};