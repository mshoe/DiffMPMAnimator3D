#include "pch.h"
#include "PointCloud.h"

std::vector<Vec3> PointCloud::GetPointPositions() const
{
	std::vector<Vec3> ret;
	ret.resize(points.size());

	for (size_t p = 0; p < points.size(); p++)
	{
		ret[p] = points[p].x;
	}

	return ret;
}