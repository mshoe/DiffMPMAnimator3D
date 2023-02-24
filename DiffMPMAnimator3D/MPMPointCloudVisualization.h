#pragma once

#include "PointCloud.h"
#include "polyscope/point_cloud.h"


void UpdatePolyscopePointCloudProperties(polyscope::PointCloud** ps_point_cloud, std::shared_ptr<DiffMPMLib3D::PointCloud> mpm_point_cloud);