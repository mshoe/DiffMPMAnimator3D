#pragma once
#include "pch.h"
#include <fstream>
namespace DiffMPMLib3D {
	struct MaterialPoint
	{
		void ResetGradients();
		void WriteEntirePointToFile(std::ofstream& ofs);
		void ReadEntirePointFromFile(std::ifstream& ifs);

		// stores every intermediate value in computation graph
		Vec3 x = Vec3::Zero();
		Vec3 v = Vec3::Zero();
		Mat3 F = Mat3::Identity();
		Mat3 C = Mat3::Zero();
		Mat3 P = Mat3::Zero();

		double m = 0.0;
		double vol = 0.0;

		double lam = 0.0;
		double mu = 0.0;

		// control variable
		Mat3 dFc = Mat3::Zero(); // How much to change deformation gradient in this timestep




		// gradients
		Vec3 dLdx = Vec3::Zero();
		Vec3 dLdv = Vec3::Zero();
		Mat3 dLdF = Mat3::Zero();
		Mat3 dLdC = Mat3::Zero();
		Mat3 dLdP = Mat3::Zero();

		Vec3 dLdv_next = Vec3::Zero();
		Mat3 dLdC_next = Mat3::Zero();

		double dLdm = 0.0;
		double dLdvol = 0.0;
	};
}