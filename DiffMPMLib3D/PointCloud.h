#pragma once
#include "pch.h"

#include "MaterialPoint.h"

namespace DiffMPMLib3D {
	struct PointCloud
	{
		void ResetGradients();
		double Compute_dLdF_Norm();
		void Descend_dLdF(double alpha, double gradient_norm);
		void RemovePoint(size_t point_index);

		void WriteToOBJ(std::string obj_path); // just writes vertex positions
		void WriteMassVelocityDefgradsToFile(std::string file_path);
		void WriteEntirePointCloudToFile(std::string file_path);
		void ReadEntirePointCloudFromFile(std::string file_path);

		std::vector<Vec3> GetPointPositions() const;
		std::vector<MaterialPoint> points;
	};
}