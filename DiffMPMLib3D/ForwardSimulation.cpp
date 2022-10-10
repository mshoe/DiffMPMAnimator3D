#include "pch.h"
#include "ForwardSimulation.h"
#include "Elasticity.h"
#include "Interpolation.h"


void ForwardTimeStep(PointCloud& next_point_cloud, PointCloud& curr_point_cloud, Grid& grid, double dt, double drag, Vec3 f_ext)
{
	P_op_1(curr_point_cloud);
	P2G(curr_point_cloud, grid, dt, drag);
	G_op(grid, dt, f_ext);
	G2P(next_point_cloud, curr_point_cloud, grid);
	P_op_2(next_point_cloud, curr_point_cloud, dt);
}

void P_op_1(PointCloud& curr_point_cloud)
{
	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		MaterialPoint& mp = curr_point_cloud.points[p];

		mp.P = PK_FixedCorotatedElasticity(mp.F + mp.dFc, mp.lam, mp.mu);
	}
}

void P2G(const PointCloud& curr_point_cloud, Grid& grid, double dt, double drag)
{
	double dx = grid.dx;

	// 1. Reset grid
	for (int i = 0; i < grid.dim_x; i++) {
		for (int j = 0; j < grid.dim_y; j++) {
			for (int k = 0; k < grid.dim_z; k++) {
				grid.nodes[i][j][k].m = 0.0;
				grid.nodes[i][j][k].v.setZero();
				grid.nodes[i][j][k].p.setZero();
			}
		}
	}


	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		const MaterialPoint& mp = curr_point_cloud.points[p];

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

			if (false/*isnan(node.p[0]) || isnan(node.p[1]) || isnan(node.p[2])*/) {
				//std::cout << "nan found" << std::endl;
				std::cout << "p: " << p << std::endl;
				std::cout << "xg: " << xg.transpose() << std::endl;
				std::cout << "xp: " << xp.transpose() << std::endl;
				std::cout << "dgp: " << dgp.transpose() << std::endl;
				std::cout << "wgp: " << wgp << std::endl;
				std::cout << "node.p: " << node.p.transpose() << std::endl;
			}
		}
	}
}

void G_op(Grid& grid, double dt, Vec3 f_ext)
{
	for (size_t i = 0; i < grid.dim_x; i++) {
		for (size_t j = 0; j < grid.dim_y; j++) {
			for (size_t k = 0; k < grid.dim_z; k++) {
				GridNode& node = grid.nodes[i][j][k];

				if (node.m != 0.0) {
					node.v = node.p / node.m + dt * f_ext;
				}
			}
		}
	}
}


void G2P(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, Grid& grid)
{
	double dx = grid.dx;
	for (size_t p = 0; p < curr_point_cloud.points.size(); p++) {
		const MaterialPoint& curr_mp = curr_point_cloud.points[p];
		MaterialPoint& next_mp = next_point_cloud.points[p];

		//next_mp = curr_mp;
		next_mp.v.setZero();
		next_mp.C.setZero();
		

		auto nodes = grid.QueryPoint_CubicBSpline(curr_mp.x);

		for (size_t i = 0; i < nodes.size(); i++) {
			const GridNode& node = nodes[i];

			Vec3 xg = node.x;
			Vec3 xp = curr_mp.x;
			Vec3 dgp = xg - xp;
			double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

			// APIC
			next_mp.v += wgp * node.v;
			next_mp.C += 3.0 / (dx * dx) * wgp * node.v * dgp.transpose();
		}
	}
}

void P_op_2(PointCloud& next_point_cloud, const PointCloud& curr_point_cloud, double dt)
{
	for (size_t p = 0; p < next_point_cloud.points.size(); p++) {
		const MaterialPoint& curr_mp = curr_point_cloud.points[p];
		MaterialPoint& next_mp = next_point_cloud.points[p];

		next_mp.F = (Mat3::Identity() + dt * next_mp.C) * (curr_mp.F + curr_mp.dFc);
		next_mp.x = curr_mp.x + dt * next_mp.v;
	}
}

void CalculatePointCloudVolumes(PointCloud& curr_point_cloud, Grid& grid)
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