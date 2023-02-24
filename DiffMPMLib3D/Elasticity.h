#pragma once
#include "pch.h"
#include "Tensor3x3x3x3.h"

namespace DiffMPMLib3D {
	Tensor3x3x3x3 d_JFit_dF_FD(const Mat3& F);

	double FixedCorotatedElasticity(const Mat3& F, double lam, double mu);
	Mat3 PK_FixedCorotatedElasticity(const Mat3& F, double lam, double mu);

	Mat3 d2_FCE_psi_dF2_mult_by_dF(const Mat3& F, double lam, double mu, const Mat3& dF);




	Tensor3x3x3x3 d2_FCE_psi_dF2_FD(const Mat3& F, double lam, double mu);

	Tensor3x3x3x3 d2_FCE_psi_dF2_mult_trick(const Mat3& F, double lam, double mu);

	Tensor3x3x3x3 d_JFit_dF(double J, const Mat3& Fit);


	void CalculateLameParameters(double young_mod, double poisson, double& lam, double& mu);
}