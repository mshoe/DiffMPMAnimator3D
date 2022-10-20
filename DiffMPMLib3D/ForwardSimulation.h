#pragma once

#include "PointCloud.h"
#include "Grid.h"


namespace SingleThreadMPM {

	void SingleParticle_op_1(MaterialPoint& mp);
	void SingleParticle_to_grid(const MaterialPoint& mp, Grid& grid, double dt, double drag);
	void SingleNode_op(GridNode& node, double dt, Vec3 f_ext);
	void Grid_to_SingleParticle(MaterialPoint& next_timestep_mp, const MaterialPoint& curr_timestep_mp, const Grid& grid);
	void SingleParticle_op_2(MaterialPoint& next_timestep_mp, const MaterialPoint& curr_timestep_mp, double dt);

	void ForwardTimeStep(PointCloud& next_point_cloud, PointCloud& curr_point_cloud, Grid& grid, double dt, double drag, Vec3 f_ext);

	void P_op_1(PointCloud& curr_point_cloud);

	void G_Reset(Grid& grid);

	void P2G(const PointCloud& curr_point_cloud, Grid& grid, double dt, double drag);

	void G_op(Grid& grid, double dt, Vec3 f_ext);

	void G2P(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, Grid& grid);

	void P_op_2(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, double dt);




	// UTILITY
	void CalculatePointCloudVolumes(PointCloud& curr_point_cloud, Grid& grid);

	void P2G_Mass(const std::vector<Vec3> points, Grid& grid, double mp_m);

	void G2G_Mass(Grid& grid_1, Grid& grid_2);

}
