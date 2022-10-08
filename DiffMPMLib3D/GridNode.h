#pragma once
#include "pch.h"


struct GridNode
{
	void ResetGradients();

	Vec3 x = Vec3::Zero(); // will never change during mpm simulation, stored for convenience i guess
	Vec3 v = Vec3::Zero();
	Vec3 p = Vec3::Zero();
	double m = 0.0;

	// derivatives
	Vec3 dLdx = Vec3::Zero();
	Vec3 dLdv = Vec3::Zero();
	Vec3 dLdp = Vec3::Zero();
	double dLdm = 0.0;
};