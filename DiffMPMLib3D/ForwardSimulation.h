#pragma once

#include "PointCloud.h"
#include "Grid.h"

namespace DiffMPMLib3D {
	namespace SingleThreadMPM {

//#define P2G_CUBIC_INTERPOLATION_LOOP(CODE) \
//Vec3 relative_point = mp.x - grid.min_point;\
//int bot_left_index[3];\
//for (size_t i = 0; i < 3; i++) {\
//	bot_left_index[i] = (int)std::floor(relative_point[i] / grid.dx) - 1;\
//}\
//for (int i = 0; i <= 3; i++) {\
//	for (int j = 0; j <= 3; j++) {\
//		for (int k = 0; k <= 3; k++) {\
//			int index[3] = {\
//				bot_left_index[0] + i,\
//				bot_left_index[1] + j,\
//				bot_left_index[2] + k\
//			};\
//			if (0 <= index[0] && index[0] < grid.dim_x &&\
//				0 <= index[1] && index[1] < grid.dim_y &&\
//				0 <= index[2] && index[2] < grid.dim_z)\
//			{\
//				const GridNode& node = grid.nodes[index[0]][index[1]][index[2]];\
//				Vec3 xg = node.x;\
//				Vec3 xp = mp.x;\
//				Vec3 dgp = xg - xp;\
//				double wgp = CubicBSpline(dgp[0] / grid.dx) * CubicBSpline(dgp[1] / grid.dx) * CubicBSpline(dgp[2] / grid.dx);\
//				CODE\
//			}\
//		}\
//	}\
//}\

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
}