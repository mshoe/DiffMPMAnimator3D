#include "pch.h"
#include "PointCloud.h"

void PointCloud::ResetGradients()
{
	for (size_t p = 0; p < points.size(); p++) {
		points[p].ResetGradients();
	}
}

double PointCloud::Compute_dLdF_Norm()
{
	double norm = 0.0;
	for (MaterialPoint& mp : points) {
		norm += mp.dLdF.squaredNorm();
	}
	norm = sqrt(norm);
	return norm;
}

void PointCloud::Descend_dLdF(double alpha, double gradient_norm)
{
	for (MaterialPoint& mp : points) 
	{
		Mat3 dir = mp.dLdF / gradient_norm;
		mp.dFc -= alpha * dir;
	}
}

std::vector<Vec3> PointCloud::GetPointPositions() const
{
	std::vector<Vec3> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = points[p].x;
	}

	return ret;
}