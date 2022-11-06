#include "pch.h"
#include "MaterialPoint.h"

void DiffMPMLib3D::MaterialPoint::ResetGradients()
{
	dLdx.setZero();
	dLdv.setZero();
	dLdF.setZero();
	dLdC.setZero();
	dLdP.setZero();

	dLdv_next.setZero();
	dLdC_next.setZero();

	dLdm = 0.0;
	dLdvol = 0.0;
}
