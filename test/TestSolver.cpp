#include "gtest/gtest.h"
#include <math.h>

#include "Solver.h"

using namespace Maxwell;

namespace AnalyticalFunctions {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;
	std::size_t standingWaveModeX = 1, standingWaveModeY = 1;

	const double PI = atan(1.0) * 4;

	double gaussianFunction(const Solver::Position& pos)
	{
		mfem::Vector normalizedPos(2);
		for (size_t i = 0; i < 2; i++) {
			double center = (meshBoundingBoxMin[i] + meshBoundingBoxMax[i]) * 0.5;
			normalizedPos[i] = 2 * (pos[i] - center) / (meshBoundingBoxMax[i] - meshBoundingBoxMin[i]);
		}

		return exp(-10. * (pow(normalizedPos[0], 2) + pow(normalizedPos[1], 2)));
	}

	double standingWaveFunction(const Solver::Position& pos)
	{
		mfem::Vector normalizedPos(2);
		mfem::Vector L(2);
		for (size_t i = 0; i < 2; i++) {
			double center = (meshBoundingBoxMin[i] + meshBoundingBoxMax[i]) * 0.5;
			L[i] = meshBoundingBoxMax[i] - meshBoundingBoxMin[i];
			normalizedPos[i] = (pos[i] - meshBoundingBoxMin[i]) / L[i];
		}

		return sin(normalizedPos[0]  * PI * double(standingWaveModeX)) *
			sin(normalizedPos[1] * PI * double(standingWaveModeY));
	}
}


class TestSolver : public ::testing::Test {
protected:
	
	std::vector<int> mapQuadElementTopLeftVertex(const mfem::Mesh& mesh) 
	{
		std::vector<int> res;
		for (int i = 0; i < mesh.GetNE(); i++) {
			mfem::Array<int> meshArrayElement;
			mesh.GetElementVertices(i, meshArrayElement);
			res.push_back(meshArrayElement[0]);
		}
		return res;
	}
};

TEST_F(TestSolver, checkRun)
{
	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	
	Solver::Options opts;
	opts.order = 1;
	opts.t_final = 1;
	opts.vis_steps = 1;
	AnalyticalFunctions::standingWaveModeX = 1;
	AnalyticalFunctions::standingWaveModeY = 1;
	
	Solver solver(opts, mesh);
	solver.getMesh().GetBoundingBox(
		AnalyticalFunctions::meshBoundingBoxMin, 
		AnalyticalFunctions::meshBoundingBoxMax);

	
	solver.setInitialElectricField(AnalyticalFunctions::standingWaveFunction);
	solver.run();
}
TEST_F(TestSolver, checkMeshDimensions) 
{

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	EXPECT_EQ(nx*ny, mesh.GetNE());

}

TEST_F(TestSolver, checkMeshElementVertices) 
{

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	std::vector<int> firstElementVerticesVector = { 0, 1, nx+2, nx+1 };
	std::vector<int> lastElementVerticesVector = { nx-1, nx, nx*2+1, nx*2 };
	mfem::Array<int> meshArrayFirstElement;
	mfem::Array<int> meshArrayLastElement;

	mesh.GetElementVertices(0, meshArrayFirstElement);
	mesh.GetElementVertices(nx*ny-1, meshArrayLastElement);

	std::vector<int> vectorFirstElement(meshArrayFirstElement.begin(), meshArrayFirstElement.end());
	std::vector<int> vectorLastElement(meshArrayLastElement.begin(), meshArrayLastElement.end());

	EXPECT_EQ(firstElementVerticesVector, vectorFirstElement);
	EXPECT_EQ(lastElementVerticesVector, vectorLastElement);

}

TEST_F(TestSolver, mapMeshElementAndVertex) 
{

	int nx = 5; int ny = 5; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	std::vector<int> mapped = mapQuadElementTopLeftVertex(mesh);

	EXPECT_EQ(0, mapped[0]);
	EXPECT_EQ(nx*ny-1, mapped.size()-1);
	EXPECT_EQ(nx-1, mapped[mapped.size()-1]);
}

TEST_F(TestSolver, checkMeshInvariance) 
{

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	Solver solver(Solver::Options(), mesh);

	std::vector<int> meshMap = mapQuadElementTopLeftVertex(mesh);
	std::vector<int> solverMeshMap = mapQuadElementTopLeftVertex(solver.getMesh());

	ASSERT_EQ(mesh.Dimension(), solver.getMesh().Dimension());
	EXPECT_EQ(meshMap[0], solverMeshMap[0]);
	EXPECT_EQ(meshMap.size(), solverMeshMap.size());
}


