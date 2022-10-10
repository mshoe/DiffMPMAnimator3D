#include "pch.h"
#include "PointCloud.h"
#include <fstream>

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

void PointCloud::WriteToOBJ(std::string obj_path)
{
	std::ofstream ofs;
	ofs.open(obj_path);

	if (ofs.good()) {
		for (size_t v = 0; v < points.size(); v++) {
			ofs << "v " << points[v].x[0] << " " << points[v].x[1] << " " << points[v].x[2] << "\n";
		}
	}
	ofs.close();
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