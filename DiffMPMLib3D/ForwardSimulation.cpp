#include "pch.h"
#include "ForwardSimulation.h"
#include "Elasticity.h"
#include "Interpolation.h"


void SingleThreadMPM::SingleParticle_op_1(MaterialPoint& mp)
{
	mp.P = PK_FixedCorotatedElasticity(mp.F + mp.dFc, mp.lam, mp.mu);
}

void SingleThreadMPM::SingleParticle_to_grid(const MaterialPoint& mp, Grid& grid, double dt, double drag)
{
	double dx = grid.dx;
	auto nodes = grid.QueryPoint_CubicBSpline(mp.x);

	for (size_t i = 0; i < nodes.size(); i++) {
		GridNode& node = nodes[i];

		Vec3 xg = node.x;
		Vec3 xp = mp.x;
		Vec3 dgp = xg - xp;
		double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

		// MLS-MPM APIC
		node.m += wgp * mp.m;
		node.p += wgp * (mp.m * mp.v * (1.0 - dt * drag) + (-3.0 / (dx * dx) * dt * mp.vol * mp.P * (mp.F + mp.dFc).transpose() + mp.m * mp.C) * dgp);
	}
}

void SingleThreadMPM::SingleNode_op(GridNode& node, double dt, Vec3 f_ext)
{
	if (node.m != 0.0) {
		node.v = node.p / node.m + dt * f_ext;
	}
}

void SingleThreadMPM::Grid_to_SingleParticle(MaterialPoint& next_timestep_mp, const MaterialPoint& curr_timestep_mp, const Grid& grid)
{
	double dx = grid.dx;
	next_timestep_mp.v.setZero();
	next_timestep_mp.C.setZero();


	auto nodes = grid.QueryPointConst_CubicBSpline(curr_timestep_mp.x);

	for (size_t i = 0; i < nodes.size(); i++) {
		const GridNode& node = nodes[i];

		Vec3 xg = node.x;
		Vec3 xp = curr_timestep_mp.x;
		Vec3 dgp = xg - xp;
		double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

		// APIC
		next_timestep_mp.v += wgp * node.v;
		next_timestep_mp.C += 3.0 / (dx * dx) * wgp * node.v * dgp.transpose();
	}
}

void SingleThreadMPM::SingleParticle_op_2(MaterialPoint& next_timestep_mp, const MaterialPoint& curr_timestep_mp, double dt)
{
	next_timestep_mp.F = (Mat3::Identity() + dt * next_timestep_mp.C) * (curr_timestep_mp.F + curr_timestep_mp.dFc);
	next_timestep_mp.x = curr_timestep_mp.x + dt * next_timestep_mp.v;
}

void SingleThreadMPM::ForwardTimeStep(PointCloud& next_point_cloud, PointCloud& curr_point_cloud, Grid& grid, double dt, double drag, Vec3 f_ext)
{
	P_op_1(curr_point_cloud);
	G_Reset(grid);
	P2G(curr_point_cloud, grid, dt, drag);
	G_op(grid, dt, f_ext);
	G2P(next_point_cloud, curr_point_cloud, grid);
	P_op_2(next_point_cloud, curr_point_cloud, dt);
}

void SingleThreadMPM::P_op_1(PointCloud& curr_point_cloud)
{
	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		MaterialPoint& mp = curr_point_cloud.points[p];

		SingleParticle_op_1(mp);
	}
}

void SingleThreadMPM::G_Reset(Grid& grid)
{
	for (int i = 0; i < grid.dim_x; i++) {
		for (int j = 0; j < grid.dim_y; j++) {
			for (int k = 0; k < grid.dim_z; k++) {
				grid.nodes[i][j][k].m = 0.0;
				grid.nodes[i][j][k].v.setZero();
				grid.nodes[i][j][k].p.setZero();
			}
		}
	}
}

void SingleThreadMPM::P2G(const PointCloud& curr_point_cloud, Grid& grid, double dt, double drag)
{
	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		const MaterialPoint& mp = curr_point_cloud.points[p];

		SingleParticle_to_grid(mp, grid, dt, drag);
	}
}

void SingleThreadMPM::G_op(Grid& grid, double dt, Vec3 f_ext)
{
	for (size_t i = 0; i < grid.dim_x; i++) {
		for (size_t j = 0; j < grid.dim_y; j++) {
			for (size_t k = 0; k < grid.dim_z; k++) {
				GridNode& node = grid.nodes[i][j][k];

				SingleNode_op(node, dt, f_ext);
			}
		}
	}
}


void SingleThreadMPM::G2P(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, Grid& grid)
{
	double dx = grid.dx;
	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		const MaterialPoint& curr_mp = curr_point_cloud.points[p];
		MaterialPoint& next_mp = next_point_cloud.points[p];

		Grid_to_SingleParticle(next_mp, curr_mp, grid);

		
	}
}

void SingleThreadMPM::P_op_2(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, double dt)
{
	for (size_t p = 0; p < next_point_cloud.points.size(); p++) {
		const MaterialPoint& curr_mp = curr_point_cloud.points[p];
		MaterialPoint& next_mp = next_point_cloud.points[p];

		SingleParticle_op_2(next_mp, curr_mp, dt);
	}
}

void SingleThreadMPM::CalculatePointCloudVolumes(PointCloud& curr_point_cloud, Grid& grid)
{
	double dx = grid.dx;
	P2G(curr_point_cloud, grid, 0.0, 0.0);

	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		MaterialPoint& mp = curr_point_cloud.points[p];

		auto nodes = grid.QueryPoint_CubicBSpline(mp.x);

		double mass_from_grid = 0.0;
		for (size_t i = 0; i < nodes.size(); i++) {
			GridNode& node = nodes[i];

			Vec3 xg = node.x;
			Vec3 xp = mp.x;
			Vec3 dgp = xg - xp;
			double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

			mass_from_grid += wgp * node.m;
		}

		mp.vol = (mass_from_grid != 0.0) ? mp.m / mass_from_grid : 0.0;
	}
}

void SingleThreadMPM::P2G_Mass(const std::vector<Vec3> points, Grid& grid, double mp_m)
{
	double dx = grid.dx;

	// 1. Reset grid
	for (int i = 0; i < grid.dim_x; i++) {
		for (int j = 0; j < grid.dim_y; j++) {
			for (int k = 0; k < grid.dim_z; k++) {
				grid.nodes[i][j][k].m = 0.0;
			}
		}
	}


	for (size_t p = 0; p < points.size(); p++) {
		const Vec3& mp_x = points[p];

		auto nodes = grid.QueryPoint_CubicBSpline(mp_x);

		for (size_t i = 0; i < nodes.size(); i++) {
			GridNode& node = nodes[i];

			Vec3 xg = node.x;
			Vec3 xp = mp_x;
			Vec3 dgp = xg - xp;
			double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

			// MLS-MPM APIC
			node.m += wgp * mp_m;
		}
	}
}

void SingleThreadMPM::G2G_Mass(Grid& grid_1, Grid& grid_2)
{
	double dx = grid_1.dx;

	// 1. Reset grid 2
	for (int i = 0; i < grid_2.dim_x; i++) {
		for (int j = 0; j < grid_2.dim_y; j++) {
			for (int k = 0; k < grid_2.dim_z; k++) {
				grid_2.nodes[i][j][k].m = 0.0;
			}
		}
	}

	// 2. Loop through grid_2 nodes as if they were particles, get interpolated mass from grid_1
	for (int i = 0; i < grid_2.dim_x; i++) {
		for (int j = 0; j < grid_2.dim_y; j++) {
			for (int k = 0; k < grid_2.dim_z; k++) {
				double& mp_m = grid_2.nodes[i][j][k].m;
				mp_m = 0.0;
				const Vec3& mp_x = grid_2.nodes[i][j][k].x;
				

				auto nodes = grid_1.QueryPoint_CubicBSpline(mp_x);
				for (size_t i = 0; i < nodes.size(); i++) {
					const GridNode& node = nodes[i];

					Vec3 xg = node.x;
					Vec3 xp = mp_x;
					Vec3 dgp = xg - xp;
					double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

					// MLS-MPM APIC
					mp_m += wgp * node.m;
				}
			}
		}
	}
}