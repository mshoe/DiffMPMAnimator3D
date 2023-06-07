#pragma once
#include "pch.h"

namespace DiffMPMLib3D
{
	void SphereUnionMarchingCubesSurfaceFromPointCloud(const std::vector<Vec3>& _points,
		double radius, double grid_dx, double iso_mass, 
		int blur_iterations,
		Vec3 grid_min_point, Vec3 grid_max_point,
		Eigen::MatrixXd& mcV, // mesh vertices
		Eigen::MatrixXi& mcF // mesh faces
	);
}