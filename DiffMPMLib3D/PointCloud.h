#pragma once
#include "pch.h"

#include "MaterialPoint.h"

struct PointCloud 
{

	std::vector<Vec3> GetPointPositions() const;
	std::vector<MaterialPoint> points;
};