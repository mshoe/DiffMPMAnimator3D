#pragma once
#include "pch.h"

#include "MaterialPoint.h"

struct PointCloud 
{
	void ResetGradients();
	double Compute_dLdF_Norm();

	std::vector<Vec3> GetPointPositions() const;
	std::vector<MaterialPoint> points;
};