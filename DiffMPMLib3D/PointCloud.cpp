#include "pch.h"
#include "PointCloud.h"
#include "Elasticity.h"
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
		ofs.precision(17);
		ofs << points.size() << std::endl;
		for (size_t v = 0; v < points.size(); v++) {
			points[v].WriteEntirePointToFile(ofs);
		}
	}
	ofs.close();
}

bool DiffMPMLib3D::PointCloud::ReadEntirePointCloudFromFile(std::string file_path)
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

		ifs.close();
		return true;
	}
	else {
		std::cout << "couldn't read file: " << file_path << std::endl;
		return false;
	}
}

void DiffMPMLib3D::PointCloud::WriteEntirePointCloudToBinaryFile(std::string file_path)
{
	std::ofstream ofs;
	ofs.open(file_path, std::ios::binary);

	if (ofs.good()) {
		cereal::BinaryOutputArchive oarchive(ofs); // Create an output archive

		oarchive(*this);
	}
	ofs.close();
}

bool DiffMPMLib3D::PointCloud::ReadEntirePointCloudFromBinaryFile(std::string file_path)
{
	std::ifstream ifs;
	ifs.open(file_path, std::ios::binary);
	if (ifs.good()) {
		cereal::BinaryInputArchive iarchive(ifs); // Create an input archive

		iarchive(*this);
	}
	else {
		std::cout << "couldn't read file: " << file_path << std::endl;
		return false;
	}
	ifs.close();
	return true;
}

bool DiffMPMLib3D::PointCloud::IsEqualToOtherPointCloud(const PointCloud& other_pc)
{
	for (size_t i = 0; i < points.size(); i++) {
		const MaterialPoint& mp = points[i];
		const MaterialPoint& other_mp = other_pc.points[i];
		if (!mp.IsEqualToOtherPoint(other_mp)) {
			
			std::cout << "Not equal at point " << i << std::endl;
			std::cout << "*** original MP ***" << std::endl;
			mp.PrintMP();
			std::cout << "*** other MP ***" << std::endl;
			other_mp.PrintMP();
			return false;
		}
	}

	std::cout << "Point clouds are equal" << std::endl;
	return true;
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

std::vector<double> DiffMPMLib3D::PointCloud::GetPointMasses() const
{
	std::vector<double> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = points[p].m;
	}

	return ret;
}

std::vector<double> DiffMPMLib3D::PointCloud::GetPointVolumes() const
{
	std::vector<double> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = points[p].vol;
	}

	return ret;
}

std::vector<DiffMPMLib3D::Mat3> DiffMPMLib3D::PointCloud::GetPointDefGrads() const
{
	std::vector<Mat3> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = points[p].F;
	}

	return ret;
}

std::vector<double> DiffMPMLib3D::PointCloud::GetPointDeterminants() const
{
	std::vector<double> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = (points[p].F + points[p].dFc).determinant();
	}

	return ret;
}

std::vector<double> DiffMPMLib3D::PointCloud::GetPointElasticEnergies() const
{
	std::vector<double> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = FixedCorotatedElasticity(points[p].F + points[p].dFc, points[p].lam, points[p].mu);
	}

	return ret;
}
