#include "pch.h"
#include "GridNode.h"

void GridNode::ResetGradients()
{
	dLdx.setZero();
	dLdv.setZero();
	dLdp.setZero();
	dLdm = 0.0;
}
