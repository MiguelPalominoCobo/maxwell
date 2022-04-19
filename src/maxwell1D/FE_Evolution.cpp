#include "FE_Evolution.h"

namespace maxwell1D {

FE_Evolution::FE_Evolution(FiniteElementSpace* fes, Options options) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	MInv_(buildInverseMassMatrix()),
<<<<<<< HEAD
	KxE_(buildDerivativeAndFluxOperator(X, Electric)),
	KxH_(buildDerivativeAndFluxOperator(X, Magnetic)),
	SxE_(buildDerivativeOperator(X, Electric)),
	SxH_(buildDerivativeOperator(X, Magnetic)),
	FxE_(buildFluxOperator(X, Electric)),
	FxH_(buildFluxOperator(X, Magnetic))
{}
=======
	K_(buildDerivativeOperator()),
	FE_(buildFluxOperators(FieldType::Electric)),
	FH_(buildFluxOperators(FieldType::Magnetic))
{
}
>>>>>>> 5b5d6a8da1ad6df994b51eb5266262e8bdfddd55

std::unique_ptr<BilinearForm> FE_Evolution::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	
	MInv->Assemble();
	MInv->Finalize();
	
	return MInv;
}

std::unique_ptr<BilinearForm> FE_Evolution::buildDerivativeOperator() const
{
	std::size_t d = 0;
	ConstantCoefficient coeff(1.0);

	auto K = std::make_unique<BilinearForm>(fes_);
	K->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(coeff, d)
		)
	);

	K->Assemble();
	K->Finalize();
	
	return K;
}

FE_Evolution::FluxOperators FE_Evolution::buildFluxOperators(const FieldType& f) const
{
	FluxOperators res = std::make_pair(
		std::make_unique<BilinearForm>(fes_),
		std::make_unique<BilinearForm>(fes_)
	);

	VectorConstantCoefficient n(Vector({ 1.0 }));
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res.first->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryFluxCoefficient(f);
		res.first->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = interiorAltFluxCoefficient();
		res.second->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryAltFluxCoefficient(f);
		res.second->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}

	res.first->Assemble();
	res.first->Finalize();
	res.second->Assemble();
	res.second->Finalize();

	return res;
}

FE_Evolution::FluxCoefficient FE_Evolution::interiorFluxCoefficient() const
{
	return FluxCoefficient{1.0, 0.0};
}

FE_Evolution::FluxCoefficient FE_Evolution::interiorAltFluxCoefficient() const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{0.0, 0.0};
	case FluxType::Upwind:
		return FluxCoefficient{0.0, -0.5}; 
	}
}

FE_Evolution::FluxCoefficient FE_Evolution::boundaryFluxCoefficient(const FieldType& f) const
{
	switch (opts_.bdrCond) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::Electric:
			return FluxCoefficient{ 0.0, 0.0 };
		case FieldType::Magnetic:
			return FluxCoefficient{ 2.0, 0.0 };
		}
	case BdrCond::SMA:
		switch (f) {
		case FieldType::Electric:
			return FluxCoefficient{ 1.0,0.0 };
		case FieldType::Magnetic:
			return FluxCoefficient{ 1.0,0.0 };
		}
	}
}

FE_Evolution::FluxCoefficient FE_Evolution::boundaryAltFluxCoefficient(const FieldType& f) const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0,0.0 };
	case FluxType::Upwind:
		switch (opts_.bdrCond) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::Electric:
				return FluxCoefficient{ 0.0, 0.0 }; // TODO
			case FieldType::Magnetic:
				return FluxCoefficient{ 0.0, 0.0 }; // TODO
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::Electric:
				return FluxCoefficient{ 0.0,0.0 };
			case FieldType::Magnetic:
				return FluxCoefficient{ 0.0,0.0 };
			}
		}
	}
}


void FE_Evolution::constructBilinearForms()
{
	MInv_ = buildInverseMassMatrix();
	K_ = buildDerivativeOperator();
	FE_ = buildFluxOperators(FieldType::Electric);
	FH_ = buildFluxOperators(FieldType::Magnetic);
}


std::unique_ptr<BilinearForm> FE_Evolution::buildDerivativeOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto S = std::make_unique<BilinearForm>(fes_);

	ConstantCoefficient one(1.0);

	S->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, d)));
	S->Assemble();
	S->Finalize();

	return S;
}


std::unique_ptr<BilinearForm> FE_Evolution::buildFluxOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto F = std::make_unique<BilinearForm>(fes_);

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	double alpha;
	double beta;

	if (ft == Electric)
	{
		alpha = -1.0;
		beta = 0.0;
		F->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		F->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}
	else
	{
		alpha = -1.0;
		beta = 0.0;
		F->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		F->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}

	F->Assemble();
	F->Finalize();

	return F;
}
void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(), fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	Vector auxRHS(MInv_->Height());

	// Update E. dE/dt = M^{-1} * (K * H - FE * {H} + altFE * [E])).
	//auxRHS = 0.0;
	K_->Mult(hOld, auxRHS);
	FE_.first->AddMult(hOld, auxRHS, -1.0);
	FE_.second->AddMult(eOld, auxRHS);
	MInv_->Mult(auxRHS, eNew);

	// Update H. dH/dt = M^{-1} * (K * E - FH * {E} + altFH * [H])).
	//auxRHS = 0.0;
	K_->Mult(eOld, auxRHS);
	FH_.first->AddMult(eOld, auxRHS, -1.0);
	FH_.second->AddMult(hOld, auxRHS);
	MInv_->Mult(auxRHS, hNew);

}

}

