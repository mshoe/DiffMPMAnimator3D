#include "pch.h"
#include "PointCloud.h"
#include <fstream>

void DiffMPMLib3D::PointCloud::ResetGradients()
{
	for (size_t p = 0; p < points.size(); p++) {
		points[p].ResetGradients();
	}
}

double DiffMPMLib3D::PointCloud::Compute_dLdF_Norm()
{
	double norm = 0.0;
	for (MaterialPoint& mp : points) {
		norm += mp.dLdF.squaredNorm();
	}
	norm = sqrt(norm);
	return norm;
}

void DiffMPMLib3D::PointCloud::Descend_dLdF(double alpha, double gradient_norm)
{
	for (MaterialPoint& mp : points) 
	{
		Mat3 dir = mp.dLdF / gradient_norm;
		mp.dFc -= alpha * dir;
	}
}

void DiffMPMLib3D::PointCloud::RemovePoint(size_t point_index)
{
	points.erase(points.begin()+point_index);
}

void DiffMPMLib3D::PointCloud::WriteToOBJ(std::string obj_path)
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

void DiffMPMLib3D::PointCloud::WriteMassVelocityDefgradsToFile(std::string file_path)
{
	std::ofstream ofs;
	ofs.open(file_path);

	if (ofs.good()) {
		for (size_t v = 0; v < points.size(); v++) {
			ofs << "x " << points[v].x[0] << " " << points[v].x[1] << " " << points[v].x[2] << "\n";
			ofs << "v " << points[v].v[0] << " " << points[v].v[1] << " " << points[v].v[2] << "\n";


			// I don't plan on doing any visualizations that actually use the 3x3 deformation gradients,
			// just maybe their elastic energy and their magnitudes?

			ofs << "F";
			for (size_t i = 0; i < 9; i++) {
				ofs << " " << points[v].F(i);
			}
			ofs << "\n";
			ofs << "dFc";
			for (size_t i = 0; i < 9; i++) {
				ofs << " " << points[v].dFc(i);
			}
			ofs << "\n";
			ofs << "dLdF";
			for (size_t i = 0; i < 9; i++) {
				ofs << " " << points[v].dLdF(i);
			}
			ofs << "\n";
		}
	}
	ofs.close();
}

void DiffMPMLib3D::PointCloud::WriteEntirePointCloudToFile(std::string file_path)
{
	// not sure how much I want to invest into this codebase in the future, so just doing the easy thing here
	std::ofstream ofs;
	ofs.open(file_path);

	if (ofs.good()) {
		// first line: number of points
		ofs << points.size() << std::endl;
		for (size_t v = 0; v < points.size(); v++) {
			points[v].WriteEntirePointToFile(ofs);
		}
	}
	ofs.close();
}

void DiffMPMLib3D::PointCloud::ReadEntirePointCloudFromFile(std::string file_path)
{
	std::ifstream ifs;
	ifs.open(file_path);

	if (ifs.good()) {
		// first line: number of points
		size_t num_points;
		ifs >> num_points;

		points.resize(num_points);

		for (size_t v = 0; v < num_points; v++) {
			points[v].ReadEntirePointFromFile(ifs);
		}
	}

	ifs.close();
}


std::vector<DiffMPMLib3D::Vec3> DiffMPMLib3D::PointCloud::GetPointPositions() const
{
	std::vector<Vec3> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = points[p].x;
	}

	return ret;
}