#pragma once
#include "pch.h"


namespace GeometryLoading 
{

	std::vector<Vec3> GeneratePointCloudFromWatertightTriangleMesh(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Vec3 min_point,
		Vec3 max_point,
		double sampling_dx
	);
}