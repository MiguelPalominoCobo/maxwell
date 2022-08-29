#pragma once

#include <functional>
#include "Types.h"
#include "Model.h"
//#include "Hopfion.h"

namespace maxwell {


class Source {
public:
	Source(Model& model, const FieldType& ft, const Direction& d,  const double spread, const double coeff,
		const Vector devFromCenter, const SourceType& srcType);

	double evalGaussianFunction3D(const Position& pos) const;
	double evalGaussianFunction2D(const Position& pos) const;
	double evalGaussianFunction1D(const Position& pos) const;
	FieldType getFieldType() const { return fieldType_; }
	SourceType getSourceType() const { return sourceType_; }
	Direction getDirection() const { return direction_; }
	double evalIniHopfionEX(const Position& pos) const;
	double evalIniHopfionEY(const Position& pos) const;
	double evalIniHopfionEZ(const Position& pos) const;
	double evalIniHopfionHX(const Position& pos) const;
	double evalIniHopfionHY(const Position& pos) const;
	double evalIniHopfionHZ(const Position& pos) const;


private:

	FieldType fieldType_;
	Direction direction_;
	SourceType sourceType_;
	double spread_;
	double coeff_;
	Vector minBB_, maxBB_, devFromCenter_;

	static Vector vectorAverage(const Vector& min, const Vector& max);
	const void checkInputArguments(Model& model);
	double normalizedPos3D(const Position& pos) const;
};

struct Sources {
public:

	void addSourceToVector(const Source& source) { sourceVector_.push_back(source); }
	const std::vector<Source>& getSourcesVector() const { return sourceVector_; }

private:

	std::vector<Source> sourceVector_;

};

}