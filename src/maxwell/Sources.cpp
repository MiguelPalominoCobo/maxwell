#include "Sources.h"
#include <math.h>
#include "Hopfion.h"

namespace maxwell {



Source::Source(
	Model& model,
	const FieldType& ft,
	const Direction& d, 
	const double spread, 
	const double coeff, 
	const Vector devFromCenter,
	const SourceType& srcType) :

	spread_(spread),
	coeff_(coeff),
	devFromCenter_(devFromCenter),
	fieldType_(ft),
	direction_(d),
	sourceType_(srcType) 
{
	checkInputArguments(model);

};

const void Source::checkInputArguments(Model& model)
{
	if (sourceType_ == SourceType::Gauss) {
		if (spread_ < 0.0) {
			throw std::exception("Invalid spread value.");
		}
		if (coeff_ < 0.0) {
			throw std::exception("Invalid coeff value.");
		}
		if (model.getConstMesh().Dimension() != devFromCenter_.Size()) {
			throw std::exception("Mesh dimension and devFromCenter vector size are not the same.");
		}
		model.getMesh().GetBoundingBox(minBB_, maxBB_, 0);
		for (int i = 0; i < devFromCenter_.Size(); i++) {
			if (devFromCenter_[i] < minBB_[i] || devFromCenter_[i] > maxBB_[i]) {
				throw std::exception("Deviation from center cannot be smaller than min boundary or bigger than max boundary values.");
			}
		}
	}
	if (sourceType_ == SourceType::Hopfion) {
		Hopfion hopfion(1, 1);             /////////////////Por ahora el 11, ya veremos luego
		if (model.getConstMesh().Dimension() != devFromCenter_.Size()) {
			throw std::exception("Mesh dimension and devFromCenter vector size are not the same.");
		}
		model.getMesh().GetBoundingBox(minBB_, maxBB_, 0);
		if (maxBB_[X] <= minBB_[X]) {
			throw std::exception("The maxBB is lower than minBB in X Direc");
		}
		if (maxBB_[Y] <= minBB_[Y]) {
			throw std::exception("The maxBB is lower than minBB in Y Direc");
		}
		if (maxBB_[Z] <= minBB_[Z]) {
			throw std::exception("The maxBB is lower than minBB in Z Direc");
		}
		for (int i = 0; i < devFromCenter_.Size(); i++) {
			if (abs(devFromCenter_[i]) < minBB_[i] || abs(devFromCenter_[i]) > maxBB_[i]) {
				throw std::exception("Deviation from center cannot be smaller than min boundary values.");
			}
		}
		for (int i = 0; i < devFromCenter_.Size(); i++) {
			if (abs(devFromCenter_[i]) > maxBB_[i]) {
				throw std::exception("Deviation from center cannot be bigger than max boundary values.");
			}
		}
		if (model.getConstMesh().Dimension() != 3) {
			throw std::exception("Mesh and Hopfion dimension are not the same.");
		}
	}
}

double Source::evalGaussianFunction3D(const Position& pos) const
{
	Vector center = vectorAverage(minBB_, maxBB_);
	Vector normalizedPos(pos.Size());
	normalizedPos = 0.0;
	for (int i = 0; i < normalizedPos.Size(); i++) {
		normalizedPos[i] = 2 * (pos[i] - center[i] - devFromCenter_[i]) / (maxBB_[i] - minBB_[i]);
	}
	return coeff_ * (1.0 / (pow(spread_, 2.0) * pow(2.0 * M_PI, 2.0 / 2.0))) *
		exp(-40 * (pow(normalizedPos[X], 2.0) + pow(normalizedPos[Y], 2.0) + pow(normalizedPos[Z], 2.0)) /
			(2.0 * pow(spread_, 2.0)));
}


double Source::evalGaussianFunction2D(const Position& pos) const
{
	Vector center = vectorAverage(minBB_, maxBB_);
	Vector normalizedPos(pos.Size());
	normalizedPos = 0.0;
	for (int i = 0; i < normalizedPos.Size(); i++) {
		normalizedPos[i] = 2 * (pos[i] - center[i] - devFromCenter_[i]) / (maxBB_[i] - minBB_[i]);
	}
	return coeff_ * (1.0 / (pow(spread_, 2.0) * pow(2.0 * M_PI, 2.0 / 2.0))) *
		exp(-40 * (pow(normalizedPos[X], 2.0) + pow(normalizedPos[Y], 2.0)) /
			(2.0 * pow(spread_, 2.0)));
}

double Source::evalGaussianFunction1D(const Position& pos) const
{
	double center = (minBB_[0] + maxBB_[0]) * 0.5 - devFromCenter_[0];
	double normalizedPos = 2 * (pos[0] - center) / (maxBB_[0] - minBB_[0]);
	return coeff_ * (1.0 / spread_ * sqrt(2.0 * M_PI)) *
		exp(-40 * pow(normalizedPos , 2.0) / pow(spread_, 2.0));

}
Vector Source::vectorAverage(const Vector& a, const Vector& b)
{
	Vector res = a;
	res.Add(1.0, b);
	res /= 2.0;
	return res;
}

double Source::normalizedPos3D(const Position& pos) const
{
	Vector center = vectorAverage(minBB_, maxBB_);
	Vector normalizedPos(pos.Size());
	normalizedPos = 0.0;
	for (int i = 0; i < normalizedPos.Size(); i++) {
		normalizedPos[i] = 2 * (pos[i] - center[i] - devFromCenter_[i]) / (maxBB_[i] - minBB_[i]);
	}
	return normalizedPos[X], normalizedPos[Y], normalizedPos[Z];

}

Vector Source::hopfionPos3D(const Position& pos) const
{
	Vector center = vectorAverage(minBB_, maxBB_);
	Vector hopfionPos(pos.Size());
	hopfionPos = 0.0;
	for (int i = 0; i < hopfionPos.Size(); i++) {
		hopfionPos[i] = 10 * (pos[i] - center[i])  + devFromCenter_[i];
	}
	return hopfionPos;

}

double Source::evalIniHopfionEX(const Position& pos) const
{
	Hopfion hopfion(1, 1);
	Vector hopfionPos = Source::hopfionPos3D(pos);
	return hopfion.getEX(-1.5, hopfionPos[X], hopfionPos[Y], hopfionPos[Z]);

}

double Source::evalIniHopfionEY(const Position& pos) const
{
	Hopfion hopfion(1, 1);
	Vector hopfionPos = Source::hopfionPos3D(pos);
	return hopfion.getEY(-1.5, hopfionPos[X], hopfionPos[Y], hopfionPos[Z]);
	//return hopfion.getEY(0, pos[X], pos[Y], pos[Z]);
}

double Source::evalIniHopfionEZ(const Position& pos) const
{
	Hopfion hopfion(1, 1);
	Vector hopfionPos = Source::hopfionPos3D(pos);
	return hopfion.getEZ(-1.5, hopfionPos[X], hopfionPos[Y], hopfionPos[Z]);
	//return hopfion.getEZ(0, pos[X], pos[Y], pos[Z]);
}

double Source::evalIniHopfionHX(const Position& pos) const
{
	Hopfion hopfion(1, 1);
	Vector hopfionPos = Source::hopfionPos3D(pos);
	return hopfion.getHX(-1.5, hopfionPos[X], hopfionPos[Y], hopfionPos[Z]);
	//return hopfion.getHX(0, pos[X], pos[Y], pos[Z]);
}

double Source::evalIniHopfionHY(const Position& pos) const
{
	Hopfion hopfion(1, 1);
	Vector hopfionPos = Source::hopfionPos3D(pos);
	return hopfion.getHX(-1.5, hopfionPos[X], hopfionPos[Y], hopfionPos[Z]);
	//return hopfion.getHY(0, pos[X], pos[Y], pos[Z]);
}

double Source::evalIniHopfionHZ(const Position& pos) const
{
	Hopfion hopfion(1, 1);
	Vector hopfionPos = Source::hopfionPos3D(pos);
	return hopfion.getHZ(-1.5, hopfionPos[X], hopfionPos[Y], hopfionPos[Z]);
}

}