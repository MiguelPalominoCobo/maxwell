#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"

#include "Types.h"

namespace maxwell1D {

using namespace mfem;

class FE_Evolution : public TimeDependentOperator {
public:

	struct Options {
		FluxType fluxType = FluxType::Upwind;
		BdrCond bdrCond = BdrCond::PEC;
	};

	enum class OperatorType {
		Stiffness,
		Flux,
		Penalty
	};
	static const std::size_t numberOfFieldComponents = 2;
	
	FE_Evolution(FiniteElementSpace* fes, Options options);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:
	struct FluxCoefficient {
		double alpha;
		double beta;
	};


	typedef std::pair<std::unique_ptr<BilinearForm>, std::unique_ptr<BilinearForm>> FluxOperators;
	typedef std::unique_ptr<BilinearForm> Operator;

	FiniteElementSpace* fes_;
<<<<<<< HEAD
	
	std::unique_ptr<BilinearForm> MInv_;
	std::unique_ptr<BilinearForm> KxE_;
	std::unique_ptr<BilinearForm> KxH_;
	std::unique_ptr<BilinearForm> SxE_;
	std::unique_ptr<BilinearForm> SxH_;
	std::unique_ptr<BilinearForm> FxE_;
	std::unique_ptr<BilinearForm> FxH_;
=======
	Options opts_;
>>>>>>> 5b5d6a8da1ad6df994b51eb5266262e8bdfddd55

	Operator MS_, FEE_, FEH_, FHE_, FHH_;
	
	void constructBilinearForms();

	Operator buildInverseMassMatrix() const;
	Operator buildDerivativeOperator() const;
	Operator buildFluxOperator(const FieldType&) const;
	Operator buildPenaltyOperator(const FieldType&) const;

	Operator buildMassAndStiffOperator() const;
	Operator buildMassAndFluxOperator(const FieldType&) const;
	Operator buildMassAndPenaltyOperator(const FieldType&) const;
	
	
	FluxOperators buildFluxOperators(const FieldType&) const;
	
	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorAltFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&) const;
	FluxCoefficient boundaryAltFluxCoefficient(const FieldType&) const;
>>>>>>> 5b5d6a8da1ad6df994b51eb5266262e8bdfddd55
};


}