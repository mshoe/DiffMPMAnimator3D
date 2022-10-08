#pragma once

#include "PointCloud.h"
#include "Grid.h"


void ForwardTimeStep(PointCloud& next_point_cloud, PointCloud& curr_point_cloud, Grid& grid, double dt, double drag, Vec3 f_ext);

void P_op_1(PointCloud& curr_point_cloud);

void P2G(const PointCloud& curr_point_cloud, Grid& grid, double dt, double drag);

void G_op(Grid& grid, double dt, Vec3 f_ext);

void G2P(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, Grid& grid);

void P_op_2(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, double dt);

void CalculatePointCloudVolumes(PointCloud& curr_point_cloud, Grid& grid);