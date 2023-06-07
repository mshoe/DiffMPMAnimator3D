#pragma once

#include "PointCloud.h"
#include "Grid.h"
#include "polyscope/point_cloud.h"


std::vector<double> GetPointCloudLocalMassFields(const DiffMPMLib3D::PointCloud& mpm_point_cloud, DiffMPMLib3D::Grid& grid);
void UpdatePolyscopePointCloudProperties(polyscope::PointCloud** ps_point_cloud, std::shared_ptr<const DiffMPMLib3D::PointCloud> mpm_point_cloud);
void UpdatePolyscopePointCloudMassField(polyscope::PointCloud** ps_point_cloud, std::shared_ptr<const DiffMPMLib3D::PointCloud> mpm_point_cloud,
	std::shared_ptr<DiffMPMLib3D::Grid> grid, double min_val, double max_val);