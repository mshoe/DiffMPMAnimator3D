#pragma once
#include "pch.h"


Mat3 PK_FixedCorotatedElasticity(const Mat3& F, double lam, double mu);

Mat3 d2_FCE_psi_dF2_mult_by_dF(const Mat3& F, double lam, double mu, const Mat3& dF);