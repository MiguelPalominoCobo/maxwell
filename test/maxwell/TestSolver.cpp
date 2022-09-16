#include "gtest/gtest.h"
#include <math.h>

#include "maxwell/Solver.h"
#include "../TestGlobalFunctions.h"

using namespace maxwell;

using Interval = std::pair<double, double>;

double analFinal(const Position& pos)
{
	double center = 0.5;
	double normalizedPos = 2 * (pos[0] - center);
	return 1.0 * (1.0 / 2.0 * sqrt(2.0 * M_PI)) *
		exp(-40 * pow(normalizedPos, 2.0) / pow(2.0, 2.0));
}


namespace AnalyticalFunctions1D {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;

	const double PI = atan(1.0) * 4;

	double gaussianFunction(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 2.0 * (pos[0] - center) /
			            ((meshBoundingBoxMax[0] - meshBoundingBoxMin[0]));
		
		return exp(-20. * pow(normalizedPos, 2));
	}

	double gaussianFunctionHalfWidth(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 4.0 * (pos[0] - center/2) /
			((meshBoundingBoxMax[0] - meshBoundingBoxMin[0]));

		return exp(-20. * pow(normalizedPos, 2));
	}
}

namespace HelperFunctions {

	void setAttributeIntervalMesh1D(
		const std::map<Attribute,Interval>& attToInterval,
		Mesh& mesh)
	{
		for (auto const& kv : attToInterval) {
			DenseMatrix changeAttMat(1, 2);
			changeAttMat.Elem(0, 0) = kv.second.first;
			changeAttMat.Elem(0, 1) = kv.second.second;
			Array<int> elemID;
			Array<IntegrationPoint> integPoint;
			mesh.FindPoints(changeAttMat, elemID, integPoint);

			if (elemID.begin() > elemID.end()) {
				throw std::exception("Lower Index bigger than Higher Index.");
			}
			if (elemID[1] > mesh.GetNE()) {
				throw std::exception("Declared element index bigger than Mesh Number of Elements.");
			}
			for (int i = elemID[0]; i <= elemID[1]; i++) {
				mesh.SetAttribute((int) i, (int) kv.first);
			}
		}
	}

	std::vector<int> mapQuadElementTopLeftVertex(
		const mfem::Mesh& mesh)
	{
		std::vector<int> res;
		for (int i = 0; i < mesh.GetNE(); i++) {
			mfem::Array<int> meshArrayElement;
			mesh.GetElementVertices(i, meshArrayElement);
			res.push_back(meshArrayElement[0]);
		}

		return res;
	}

	std::map<Time, FieldFrame>::const_iterator findTimeId(
		const std::map<Time, FieldFrame>& timeMap,
		const Time& timeToFind,
		const double tolerance)
	{
		for (auto it = timeMap.begin(); it != timeMap.end(); it++) {
			const Time& time = it->first;
			if (abs(time - timeToFind) < tolerance) {
				return it;
			}
		}
		return timeMap.end();
	}

}

using namespace AnalyticalFunctions1D;
using namespace HelperFunctions;

class TestMaxwellSolver : public ::testing::Test {
protected:

	Model buildOneDimOneMatModel(
		const int meshIntervals = 51, 
		const BdrCond& bdrL = BdrCond::PEC, 
		const BdrCond& bdrR = BdrCond::PEC) {

		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), AttributeToMaterial(), buildAttrToBdrMap1D(bdrL,bdrR));
	}

	Sources buildSourcesWithDefaultSource(
		Model& model, 
		const FieldType& ft = E,
		const Direction& d = X, 
		const double spread = 2.0, 
		const double coeff = 1.0, 
		const Vector dev = Vector({ 0.0 })) {

		Sources res;
		res.addSourceToVector(Source(model, ft, d,  spread, coeff, dev, SourceType::Gauss));
		return res;
	}

	Probes buildProbesWithDefaultPointsProbe(
		const FieldType& fToExtract = E, 
		const Direction& dirToExtract = X)
	{
		Probes res;
		res.vis_steps = 20;
		res.addPointsProbeToCollection(PointsProbe(fToExtract, dirToExtract,
			std::vector<std::vector<double>>{{0.0},{0.5},{1.0}}));
		return res;
	}

	Probes buildProbes2DWithDefaultPointsProbe(
		const FieldType& fToExtract = E,
		const Direction& dirToExtract = X)
	{
		Probes res;
		res.addPointsProbeToCollection(PointsProbe(fToExtract, dirToExtract,
			std::vector<std::vector<double>>{ {0.0, 0.0}, { 0.5, 0.5 }, { 1,1 }}));
		return res;
	}

	Probes buildProbesWithDefaultPointsProbeRB(
		const FieldType& fToExtract = E,
		const Direction& dirToExtract = Y)
	{
		// Resonant box 0.06 x 0.05 x 0.08
		Probes res;
		res.vis_steps = 20;
		res.addPointsProbeToCollection(PointsProbe(fToExtract, dirToExtract,
			std::vector<std::vector<double>>{ { 0.03, 0.025, 0.00 }, { 0.000, 0.025, 0.04 }, { 0.030, 0.000, 0.04 }, 
											  { 0.03, 0.025, 0.04 }, { 0.015, 0.0375, 0.06 }, { 0.045, 0.0125, 0.02 }}));

		res.addPointsProbeToCollection(PointsProbe(E, Z,
			std::vector<std::vector<double>>{ { 0.03, 0.025, 0.00 }, { 0.000, 0.025, 0.04 }, { 0.030, 0.000, 0.04 },
											  { 0.03, 0.025, 0.04 }, { 0.015, 0.0375, 0.06 }, { 0.045, 0.0125, 0.02 }}));

		res.addPointsProbeToCollection(PointsProbe(E, X,
			std::vector<std::vector<double>>{ { 0.03, 0.025, 0.00 }, { 0.000, 0.025, 0.04 }, { 0.030, 0.000, 0.04 },
											  { 0.03, 0.025, 0.04 }, { 0.015, 0.0375, 0.06 }, { 0.045, 0.0125, 0.02 }}));

		res.addPointsProbeToCollection(PointsProbe(H, Y,
			std::vector<std::vector<double>>{ { 0.03, 0.025, 0.00 }, { 0.000, 0.025, 0.04 }, { 0.030, 0.000, 0.04 },
											  { 0.03, 0.025, 0.04 }, { 0.015, 0.0375, 0.06 }, { 0.045, 0.0125, 0.02 }}));
		return res;
	}

	maxwell::Solver::Options buildDefaultSolverOpts(const double tFinal = 2.0)
	{
		maxwell::Solver::Options res;

		res.evolutionOperatorOptions = FiniteElementEvolution::Options();
		res.t_final = tFinal;

		return res;
	}

	AttributeToBoundary buildAttrToBdrMap1D(const BdrCond& bdrL, const BdrCond& bdrR)
	{
		return {
			{1, bdrL},
			{2, bdrR}
		};
	}

	AttributeToMaterial buildAttToVaccumOneMatMap1D()
	{
		return { 
			{ 1, Material(1.0, 1.0) } 
		};
	}

	double getBoundaryFieldValueAtTime(
		const PointsProbe& probe,
		const Time& timeToFind,
		const int denseMatPointByOrder)
	{
		auto itpos = HelperFunctions::findTimeId(probe.getConstFieldMovie(), timeToFind, 1e-6);
		if (itpos == probe.getConstFieldMovie().end()) {
			throw std::exception("Time value has not been found within the specified tolerance.");
		}
		auto FieldValueForTimeAtPoint = itpos->second.at(denseMatPointByOrder).at(probe.getDirection());

		return FieldValueForTimeAtPoint;
	}

	double computeEnergy(GridFunction& eX, GridFunction& eY, GridFunction& eZ, GridFunction& hX, GridFunction& hY, GridFunction& hZ) {
		return pow(eX.Norml2(), 2.0) + pow(eY.Norml2(), 2.0) + pow(eZ.Norml2(), 2.0)
				+ pow(hX.Norml2(), 2.0) + pow(hY.Norml2(), 2.0) + pow(hZ.Norml2(), 2.0);
	}
};
TEST_F(TestMaxwellSolver, oneDimensional_centered)
{	
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts, with its mesh being one-dimensional.
	The field along Ey is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that after two seconds with PEC boundary conditions, the wave evolves
	back to its initial state within the specified error.*/

	Model model = buildOneDimOneMatModel();

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(
		model, 
		Probes(),
		buildSourcesWithDefaultSource(model), 
		solverOpts);
	
	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_centered_energy)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts, with its mesh being one-dimensional.
	The field along Ey is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that after two seconds with PEC boundary conditions, the wave evolves
	back to its initial state within the specified error.*/

	Model model = buildOneDimOneMatModel();

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(
		model,
		Probes(),
		buildSourcesWithDefaultSource(model),
		solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	EXPECT_GE(pow(eOld.Norml2(),2.0) + pow(hOld.Norml2(),2.0), pow(eNew.Norml2(),2.0) + pow(hNew.Norml2(),2.0));

}
TEST_F(TestMaxwellSolver, oneDimensional_centered_PEC_EY)
{
	Model model = buildOneDimOneMatModel();

	auto probes = buildProbesWithDefaultPointsProbe(E, Y);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts());

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}

TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC_EX)
{
	Mesh mesh = Mesh::MakeCartesian1D(51, 1.0);
	Model model = Model(mesh, AttributeToMaterial(), buildAttrToBdrMap1D(BdrCond::PEC, BdrCond::PEC));

	Probes probes;
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();

	maxwell::Solver solver(
		model, 
		probes,
		buildSourcesWithDefaultSource(model),
		solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);



	DG_FECollection fec(solverOpts.order, mesh.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace fes(&mesh, &fec);
	mfem::FunctionCoefficient u0(analFinal);
	GridFunction u(&fes);
	u.ProjectCoefficient(u0);
	EXPECT_NEAR(u.Norml2(), eNew.Norml2(), 2e-3);



}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC_EY)
{
	Model model = buildOneDimOneMatModel();

	auto probes = buildProbesWithDefaultPointsProbe(E, Y);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts());

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC_EZ)
{
	Mesh mesh = Mesh::MakeCartesian1D(51, 1.0);
	Model model = Model(mesh, AttributeToMaterial(), buildAttrToBdrMap1D(BdrCond::PEC, BdrCond::PEC));

	auto probes = buildProbesWithDefaultPointsProbe(E, Z);
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Z),
		solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

	DG_FECollection fec(solverOpts.order, mesh.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace fes(&mesh, &fec);
	mfem::FunctionCoefficient u0(analFinal);
	GridFunction u(&fes);
	u.ProjectCoefficient(u0);
	EXPECT_NEAR(u.Norml2(), eNew.Norml2(), 2e-3);

	std::unique_ptr<mfem::ParaViewDataCollection> pd = std::make_unique<ParaViewDataCollection>("Maxwell", &mesh);
	pd->SetPrefixPath("ParaView");
	pd->RegisterField("EzTheo", &u);
	pd->SetLevelsOfDetail(solverOpts.order);
	pd->SetDataFormat(VTKFormat::BINARY);
	solverOpts.order > 0 ? pd->SetHighOrderOutput(true) : pd->SetHighOrderOutput(false);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();

}


TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC_HX)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::PMC, BdrCond::PMC);

	auto probes = buildProbesWithDefaultPointsProbe(H, Y);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, H, X),
		buildDefaultSolverOpts());

	GridFunction hOld = solver.getFieldInDirection(H, X);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, X);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC_HY)
{
	Mesh mesh = Mesh::MakeCartesian1D(51, 1.0);
	Model model = Model(mesh, AttributeToMaterial(), buildAttrToBdrMap1D(BdrCond::PMC, BdrCond::PMC));

	auto probes = buildProbesWithDefaultPointsProbe(H, Y);
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, H, Y),
		solverOpts
	);

	GridFunction hOld = solver.getFieldInDirection(H, Y);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Y);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

	DG_FECollection fec(solverOpts.order, mesh.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace fes(&mesh, &fec);
	mfem::FunctionCoefficient u0(analFinal);
	GridFunction u(&fes);
	u.ProjectCoefficient(u0);
	EXPECT_NEAR(u.Norml2(), hNew.Norml2(), 2e-3);

	double errorL2 = u.Norml2() - hNew.Norml2();
	EXPECT_NEAR(0.0, error, 2e-3);

	std::unique_ptr<mfem::ParaViewDataCollection> pd = std::make_unique<ParaViewDataCollection>("Maxwell", &mesh);
	pd->SetPrefixPath("ParaView");
	pd->RegisterField("HyTheo", &u);
	pd->SetLevelsOfDetail(solverOpts.order);
	pd->SetDataFormat(VTKFormat::BINARY);
	solverOpts.order > 0 ? pd->SetHighOrderOutput(true) : pd->SetHighOrderOutput(false);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC_HZ)
{
	Mesh mesh = Mesh::MakeCartesian1D(51, 1.0);
	Model model = Model(mesh, AttributeToMaterial(), buildAttrToBdrMap1D(BdrCond::PMC, BdrCond::PMC));

	auto probes = buildProbesWithDefaultPointsProbe(H, Z);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, H, Z),
		solverOpts);

	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);
																					  
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

	DG_FECollection fec(solverOpts.order, mesh.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace fes(&mesh, &fec);
	mfem::FunctionCoefficient u0(analFinal);
	GridFunction u(&fes);
	u.ProjectCoefficient(u0);
	EXPECT_NEAR(u.Norml2(), hNew.Norml2(), 2e-3);

	double errorL2 = u.Norml2() - hNew.Norml2();
	EXPECT_NEAR(0.0, error, 2e-3);

	std::unique_ptr<mfem::ParaViewDataCollection> pd = std::make_unique<ParaViewDataCollection>("Maxwell", &mesh);
	pd->SetPrefixPath("ParaView");
	pd->RegisterField("HzTheo", &u);
	pd->SetLevelsOfDetail(solverOpts.order);
	pd->SetDataFormat(VTKFormat::BINARY);
	solverOpts.order > 0 ? pd->SetHighOrderOutput(true) : pd->SetHighOrderOutput(false);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();

}

TEST_F(TestMaxwellSolver, DISABLED_oneDimensional_upwind_SMA_EX)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(E, X),
		buildSourcesWithDefaultSource(model, E, X),
		buildDefaultSolverOpts(0.2));

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, DISABLED_oneDimensional_upwind_SMA_EY)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	auto probes = buildProbesWithDefaultPointsProbe(E, Y);
	auto probeEZ = PointsProbe(E, Z, std::vector<std::vector<double>>({ {0.0},{0.5},{1.0} }));
	auto probeHY = PointsProbe(H, Y, std::vector<std::vector<double>>({ {0.0},{0.5},{1.0} }));
	probes.addPointsProbeToCollection(probeEZ);
	probes.addPointsProbeToCollection(probeHY);
	probes.addExporterProbeToCollection(ExporterProbe());
	

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts(1.0));

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 2), 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 2), 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.0, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.0, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.0, 2), 2e-3);

}
TEST_F(TestMaxwellSolver, DISABLED_oneDimensional_upwind_SMA_EZ)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	auto probes = buildProbesWithDefaultPointsProbe(E, Z);
	//probes.addExporterProbeToCollection(ExporterProbe());
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(1.0);
	solverOpts.order = 4;

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Z),
		solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_GE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_GE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2));

}

TEST_F(TestMaxwellSolver, oneDimensional_strong_flux_PEC_EY)
{
	Mesh mesh = Mesh::MakeCartesian1D(51);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	maxwell::Solver::Options opts;
	opts.evolutionOperatorOptions = FiniteElementEvolution::Options();
	opts.evolutionOperatorOptions.disForm = DisForm::Strong;
	opts.t_final = 0.5;

	Probes probes = buildProbesWithDefaultPointsProbe(E, Y);
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 5;

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		opts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}

TEST_F(TestMaxwellSolver, DISABLED_oneDimensional_weak_strong_flux_comparison)
{
	Model model = buildOneDimOneMatModel();
	
	maxwell::Solver::Options optsWeak;
	optsWeak.evolutionOperatorOptions = FiniteElementEvolution::Options();

	maxwell::Solver solverWeak(
		model,
		buildProbesWithDefaultPointsProbe(E, Y),
		buildSourcesWithDefaultSource(model, E, Y),
		optsWeak);

	maxwell::Solver::Options optsStrong;
	optsStrong.evolutionOperatorOptions = FiniteElementEvolution::Options();
	optsStrong.evolutionOperatorOptions.disForm = DisForm::Strong;

	maxwell::Solver solverStrong(
		model,
		buildProbesWithDefaultPointsProbe(E, Y),
		buildSourcesWithDefaultSource(model, E, Y),
		optsStrong);

	GridFunction eOldWk = solverWeak.getFieldInDirection(E, Y);
	GridFunction eOldSt = solverStrong.getFieldInDirection(E, Y);
	solverWeak.run();
	solverStrong.run();
	GridFunction eNewWk = solverWeak.getFieldInDirection(E, Y);
	GridFunction eNewSt = solverStrong.getFieldInDirection(E, Y);

	double errorOld = eOldWk.DistanceTo(eOldSt);
	double errorNew = eNewWk.DistanceTo(eNewSt);
	EXPECT_NEAR(0.0, errorOld, 2e-3);
	EXPECT_NEAR(0.0, errorNew, 2e-3);

}

TEST_F(TestMaxwellSolver, twoSourceWaveTravelsToTheRight_SMA)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	Probes probes;
	probes.addPointsProbeToCollection(PointsProbe(E, Y, std::vector<std::vector<double>>{ {0.5}, { 0.8 } }));
	probes.addExporterProbeToCollection(ExporterProbe());

	Sources sources;
	sources.addSourceToVector(Source(model, E, Y, 2.0, 1.0, Vector({ 0.0 }), SourceType::Gauss));
	sources.addSourceToVector(Source(model, H, Z, 2.0, 1.0, Vector({ 0.0 }), SourceType::Gauss));

	maxwell::Solver solver(
		model,
		probes,
		sources,
		buildDefaultSolverOpts(0.7));

	solver.run();

	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.3, 1),
				getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0),
				2e-3);

}
TEST_F(TestMaxwellSolver, twoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	Mesh mesh1D = Mesh::MakeCartesian1D(101);
	setAttributeIntervalMesh1D({ { 2, std::make_pair(0.76, 1.0) } }, mesh1D);

	Model model = Model(
		mesh1D, 
		{
			{1, Material(1.0, 1.0)},
			{2, Material(2.0, 1.0)}
		},
		{
			{1, BdrCond::SMA},
			{2, BdrCond::PEC}
		}
	);

	Probes probes;
	probes.addPointsProbeToCollection(PointsProbe(E, Y, std::vector<std::vector<double>>{ {0.3}, { 0.1 } }));

	Sources sources;
	sources.addSourceToVector(Source(model, E, Y, 1.0, 0.5, Vector({ 0.2 }), SourceType::Gauss));
	sources.addSourceToVector(Source(model, H, Z, 1.0, 0.5, Vector({ 0.2 }), SourceType::Gauss));

	maxwell::Solver solver(
		model,
		probes,
		sources,
		buildDefaultSolverOpts(1.5));

	auto eOld = solver.getFieldInDirection(E, Y);

	double reflectCoeff =
		(model.getAttToMat().at(2).getImpedance() - model.getAttToMat().at(1).getImpedance()) /
		(model.getAttToMat().at(2).getImpedance() + model.getAttToMat().at(1).getImpedance());

	solver.run();

	EXPECT_NEAR(eOld.Max(), 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0), 2e-3);
	EXPECT_NEAR(0.0, 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.45, 0), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.90, 0), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.10, 1), 2e-3);
	EXPECT_NEAR(0.0, 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.30, 1), 2e-3);
}
TEST_F(TestMaxwellSolver, twoDimensional_Periodic) //TODO ADD ENERGY CHECK
{
	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
	Vector periodic({ 0.0, 1.0 });
	std::vector<Vector> trans;
	trans.push_back(periodic);
	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D,mesh2D.CreatePeriodicVertexMapping(trans));
	
	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());
	
	Probes probes;
	//probes.addExporterProbeToCollection(ExporterProbe());
	//probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.2, 0.0 }), SourceType::Gauss));

	maxwell::Solver solver(model, probes, sources, buildDefaultSolverOpts(1.0));

	solver.run();

}

TEST_F(TestMaxwellSolver, DISABLED_twoDimensional_Periodic_strong) //TODO ADD ENERGY CHECK
{
	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
	Vector periodic({ 0.0, 1.0 });
	std::vector<Vector> trans;
	trans.push_back(periodic);
	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D, mesh2D.CreatePeriodicVertexMapping(trans));

	maxwell::Solver::Options opts;
	opts.evolutionOperatorOptions = FiniteElementEvolution::Options();
	opts.evolutionOperatorOptions.disForm = DisForm::Strong;

	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());

	Probes probes;
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.0, 0.0 }), SourceType::Gauss));

	maxwell::Solver solver(
		model, 
		probes, 
		sources, 
		opts);

	solver.run();

}
TEST_F(TestMaxwellSolver, twoDimensional_centered_NC_MESH) //TODO ADD ENERGY CHECK
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.
	The field along Ez is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty 
	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not 
	higher than the initial value.*/

	const char* mesh_file = "../maxwell/mesh/star-mixed.mesh";
	Mesh mesh(mesh_file);
	mesh.UniformRefinement(6);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Probes probes;
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, Z, 2.0, 1.0, Vector({ 0.0, 0.0 }), SourceType::Gauss));

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(3.0);
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.order = 6;

	maxwell::Solver solver(model, probes, sources, solverOpts);

	GridFunction eOldX = solver.getFieldInDirection(E, X);
	GridFunction eOldY = solver.getFieldInDirection(E, Y);
	GridFunction eOldZ = solver.getFieldInDirection(E, Z);
	GridFunction hOldX = solver.getFieldInDirection(H, X);
	GridFunction hOldY = solver.getFieldInDirection(H, Y);
	GridFunction hOldZ = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNewX = solver.getFieldInDirection(E, X);
	GridFunction eNewY = solver.getFieldInDirection(E, Y);
	GridFunction eNewZ = solver.getFieldInDirection(E, Z);
	GridFunction hNewX = solver.getFieldInDirection(H, X);
	GridFunction hNewY = solver.getFieldInDirection(H, Y);
	GridFunction hNewZ = solver.getFieldInDirection(H, Z);

	EXPECT_GT(eOldZ.Max(), eNewZ.Max());

	double Eie = pow(eOldX.Norml2(), 2.0) + pow(eOldY.Norml2(), 2.0) + pow(eOldZ.Norml2(), 2.0);
	double Eih = pow(hOldX.Norml2(), 2.0) + pow(hOldY.Norml2(), 2.0) + pow(hOldZ.Norml2(), 2.0);

	double Efe = pow(eNewX.Norml2(), 2.0) + pow(eNewY.Norml2(), 2.0) + pow(eNewZ.Norml2(), 2.0);
	double Efh = pow(hNewX.Norml2(), 2.0) + pow(hNewY.Norml2(), 2.0) + pow(hNewZ.Norml2(), 2.0);

	EXPECT_GE(Eie + Eih, Efe + Efh);
}
TEST_F(TestMaxwellSolver, twoDimensional_centered_AMR_MESH)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.
	The field along Ez is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
	higher than the initial value.*/

	const char* mesh_file = "../maxwell/mesh/amr-quad.mesh";
	Mesh mesh(mesh_file);
	mesh.UniformRefinement();
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Sources sources;
	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 }), SourceType::Gauss));

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(2.92);
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(model, Probes(), sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);

	EXPECT_GT(eOld.Max(), eNew.Max());
}
TEST_F(TestMaxwellSolver, checkFluxOperatorO2)
{

	Model model = buildOneDimOneMatModel(3);

	maxwell::Solver solver(
		model,
		Probes(),
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts(0.1));

	auto MSMat = convertMFEMDenseToEigen(solver.getFEEvol().get()
		->getInvMassStiffness(FieldType::E, Direction::X).get()->SpMat().ToDenseMatrix());
	auto noDirMat =  convertMFEMDenseToEigen(solver.getFEEvol().get()
		->getInvMassNoDirFlux(FieldType::E, FieldType::E).get()->SpMat().ToDenseMatrix());
	auto oneDirMat = convertMFEMDenseToEigen(solver.getFEEvol().get()
		->getInvMassOneDirFlux(FieldType::E, FieldType::H, Direction::X).get()->SpMat().ToDenseMatrix());

	std::cout << MSMat << std::endl;

}

TEST_F(TestMaxwellSolver, twoDimensional_Centered_PEC_EZ)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.

	The field along Ez is extracted before and after the solver calls its run() method
	and evolves the problem. After two seconds, this test verifies that the problem reaches 
	a new peak  in field Ez. check some points where the fields must be null and the maximum
	value in Ez is not higher than the initial value.*/

	Mesh mesh = Mesh::MakeCartesian2D(11, 11, Element::Type::QUADRILATERAL);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Sources sources;
	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 }), SourceType::Gauss));

	auto probes = buildProbes2DWithDefaultPointsProbe(E, Z);
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(model, probes, sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	GridFunction hOldX = solver.getFieldInDirection(H, X);
	GridFunction hOldY = solver.getFieldInDirection(H, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);
	GridFunction hNewX = solver.getFieldInDirection(H, X);
	GridFunction hNewY = solver.getFieldInDirection(H, Y);

	EXPECT_GT(eOld.Max(), eNew.Max());

	double Ei = pow(eOld.Norml2(), 2.0) + pow(hOldX.Norml2(), 2.0) + pow(hOldY.Norml2(), 2.0);
	double Ef = pow(eNew.Norml2(), 2.0) + pow(hNewX.Norml2(), 2.0) + pow(hNewY.Norml2(), 2.0);

	EXPECT_GE(Ei, Ef);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));
}

TEST_F(TestMaxwellSolver, threeDimensional_centered_PEC_EZ)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A PEC box is implemented in this code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 3D Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts. The fields are extracted before and 
	after the solver calls its run() method and evolves the problem.

	The field along Ez is extracted before and after the solver calls its run() method
	and evolves the problem. After two seconds, this test verifies that the problem reaches 
	a new peak  in field Ez. check some points where the fields must be null and the maximum
	value in Ez is not higher than the initial value.*/

	Mesh mesh = Mesh::MakeCartesian3D(3, 3, 3, Element::Type::HEXAHEDRON);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Probes probes;
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, Z, 2.0, 1.0, Vector({ 0.0, 0.0, 0.0 }), SourceType::Gauss));

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.order = 3;

	maxwell::Solver solver(model, probes, sources, solverOpts);

	GridFunction eOldX = solver.getFieldInDirection(E, X);
	GridFunction eOldY = solver.getFieldInDirection(E, Y);
	GridFunction eOldZ = solver.getFieldInDirection(E, Z);
	GridFunction hOldX = solver.getFieldInDirection(H, X);
	GridFunction hOldY = solver.getFieldInDirection(H, Y);
	GridFunction hOldZ = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNewX = solver.getFieldInDirection(E, X);
	GridFunction eNewY = solver.getFieldInDirection(E, Y);
	GridFunction eNewZ = solver.getFieldInDirection(E, Z);
	GridFunction hNewX = solver.getFieldInDirection(H, X);
	GridFunction hNewY = solver.getFieldInDirection(H, Y);
	GridFunction hNewZ = solver.getFieldInDirection(H, Z);


	EXPECT_GT(eOldZ.Max(), eNewZ.Max());

	double Eie = pow(eOldX.Norml2(), 2.0) + pow(eOldY.Norml2(), 2.0) + pow(eOldZ.Norml2(), 2.0);
	double Eih = pow(hOldX.Norml2(), 2.0) + pow(hOldY.Norml2(), 2.0) + pow(hOldZ.Norml2(), 2.0);

	double Efe = pow(eNewX.Norml2(), 2.0) + pow(eNewY.Norml2(), 2.0) + pow(eNewZ.Norml2(), 2.0);
	double Efh = pow(hNewX.Norml2(), 2.0) + pow(hNewY.Norml2(), 2.0) + pow(hNewZ.Norml2(), 2.0);

	EXPECT_GE(Eie + Eih, Efe + Efh);

}
TEST_F(TestMaxwellSolver, threeDimensionalResonantBox)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A rectangular resonant cavity is implemented in this code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 3D Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts. The fields are extracted before and
	after the solver calls its run() method and evolves the problem.

	The field along Ez is extracted before and after the solver calls its run() method
	and evolves the problem. After two seconds, this test verifies that the problem reaches
	a new peak  in field Ez. check that the fields must be null in certain points and the maximum
	value in Ez is not higher than the initial value.*/

	Mesh mesh = Mesh::MakeCartesian3D(6, 6, 6, Element::Type::HEXAHEDRON, 0.06, 0.05, 0.08);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Sources sources;
	sources.addSourceToVector(Source(model, E, Y, 2.0, 20.0, Vector({ 0.0, 0.0, 0.0 }), SourceType::Gauss));

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.order = 3;

	maxwell::Solver solver(model, buildProbesWithDefaultPointsProbeRB(), sources, solverOpts);

	GridFunction eOldX = solver.getFieldInDirection(E, X);
	GridFunction eOldY = solver.getFieldInDirection(E, Y);
	GridFunction eOldZ = solver.getFieldInDirection(E, Z);
	GridFunction hOldX = solver.getFieldInDirection(H, X);
	GridFunction hOldY = solver.getFieldInDirection(H, Y);
	GridFunction hOldZ = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNewX = solver.getFieldInDirection(E, X);
	GridFunction eNewY = solver.getFieldInDirection(E, Y);
	GridFunction eNewZ = solver.getFieldInDirection(E, Z);
	GridFunction hNewX = solver.getFieldInDirection(H, X);
	GridFunction hNewY = solver.getFieldInDirection(H, Y);
	GridFunction hNewZ = solver.getFieldInDirection(H, Z);


	EXPECT_GT(eOldY.Max(), eNewY.Max());

	EXPECT_GE(computeEnergy(eOldX, eOldY, eOldZ, hOldX, hOldY, hOldZ), computeEnergy(eNewX, eNewY, eNewZ, hNewX, hNewY, hNewZ));

	// Punto 1: { 0.03, 0.025, 0.00 }
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 1.5, 0), 2e-3);

	// Punto 2: { 0.000, 0.025, 0.04 }
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 0.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 0.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 0.5, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 1.5, 1), 2e-3);

	// Punto 3: { 0.030, 0.000, 0.04 }
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 1.5, 2), 2e-3);


	// Punto 4:  { 0.03, 0.025, 0.04 } = center
	EXPECT_NE(eOldY.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 3));
	EXPECT_NE(eOldY.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 3));
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 0.5, 3), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.5, 3), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 0.5, 3), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.5, 3), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 0.5, 3), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 1.5, 3), 2e-3);


	// Punto 5: { 0.015, 0.375, 0.06 }
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 0.5, 4), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.5, 4), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 0.5, 4), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.5, 4), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 0.5, 4), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 1.5, 4), 2e-3);

	// Punto 6:  { 0.045, 0.125, 0.02 }
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 0.5, 5), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.5, 5), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 0.5, 3), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.5, 5), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 0.5, 5), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(3), 1.5, 5), 2e-3);

}



//TEST_F(TestMaxwellSolver, DISABLED_twoDimensionalResonantBox)
//{
//	Mesh mesh2D = Mesh::MakeCartesian2D(21, 21, Element::Type::QUADRILATERAL);
//	std::vector<Attribute> attArrSingle = std::vector<Attribute>({ 1 });
//	Material mat11 = Material(1.0, 1.0);
//	std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
//	AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatMap(attArrSingle, matArrSimple);
//	AttributeToBoundary attToBdrVec;
//	Model model(mesh2D, attToMatVec, attToBdrVec);
//
//	double spread = 2.0;
//	double coeff = 20.0;
//	const Vector dev = Vector({ 0.0,0.0 });
//	Source EXFieldSource = Source(model, spread, coeff, dev, X, E); 
//	Sources sources;
//	sources.addSourceToVector(EXFieldSource);
//
//	Probes probes;
//	//probes.paraview = true;
//	probes.vis_steps = 100;
//
//	maxwell::Solver::Options solverOpts;
//
//	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
//	solverOpts.dt = 1e-4;
//	solverOpts.order = 1;
//
//	maxwell::Solver solver(model, probes,
//		sources, solverOpts);
//
//	solver.run();
//
//}