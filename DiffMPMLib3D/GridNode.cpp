#include "pch.h"
#include "GridNode.h"

void DiffMPMLib3D::GridNode::ResetGradients()
{
	dLdx.setZero();
	dLdv.setZero();
	dLdp.setZero();
	dLdm = 0.0;
}
