#pragma once
#include "pch.h"
#include "GridNode.h"
#include <functional>

struct Grid
{
	/* 
	* Grid stores a 3D uniform grid of nodal interpolation points.
	* The convention is that i, j, k corresponds the the positive x, y, z axes.
	*/

	Grid(int _dim_x, int _dim_y, int dim_z, double _dx, Vec3 _min_point);

	std::vector<std::reference_wrapper<GridNode>> QueryPoint_CubicBSpline(Vec3 point);

	std::vector<Vec3> GetNodePositions() const;


	int dim_x;
	int dim_y;
	int dim_z;
	double dx;
	Vec3 min_point;
	std::vector<std::vector<std::vector<GridNode>>> nodes;
};