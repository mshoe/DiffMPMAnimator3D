#pragma once
#include "pch.h"

#include "MaterialPoint.h"

#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

namespace DiffMPMLib3D {
	struct PointCloud
	{
		void ResetGradients();
		double Compute_dLdF_Norm();
		void Descend_dLdF(double alpha, double gradient_norm);
		void RemovePoint(size_t point_index);

		bool ReadFromOBJ(std::string obj_path, double point_mass);
		void WriteToOBJ(std::string obj_path); // just writes vertex positions
		void WriteMassVelocityDefgradsToFile(std::string file_path);
		/*void WriteEntirePointCloudToFile(std::string file_path);
		bool ReadEntirePointCloudFromFile(std::string file_path);*/
		void WriteEntirePointCloudToBinaryFile(std::string file_path);
		bool ReadEntirePointCloudFromBinaryFile(std::string file_path);

		bool IsEqualToOtherPointCloud(const PointCloud& other_pc);

		std::vector<Vec3> GetPointPositions() const;
		std::vector<double> GetPointMasses() const;
		std::vector<double> GetPointVolumes() const;
		std::vector<Mat3> GetPointDefGrads() const;
		std::vector<double> GetPointDeterminants() const;
		std::vector<double> GetPointElasticEnergies() const;

		std::vector<MaterialPoint> points;

		template<class Archive>
		void serialize(Archive& archive)
		{
			archive(
				CEREAL_NVP(points)
			);
		}
	};
}