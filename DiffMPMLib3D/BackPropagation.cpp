#include "pch.h"
#include "BackPropagation.h"
#include "Elasticity.h"
#include "Interpolation.h"

void DiffMPMLib3D::Back_Timestep(CompGraphLayer& layer_nplus1, CompGraphLayer& layer_n, double drag, double dt)
{
	
	const PointCloud& pc = *layer_nplus1.point_cloud;
	PointCloud& pc_prev = *layer_n.point_cloud;
	Grid& grid = *layer_n.grid;
	double dx = grid.dx;

	// Back prop P_op_2
	for (size_t p = 0; p < pc.points.size(); p++) {
		const MaterialPoint& mp = pc.points[p];
		MaterialPoint& mp_prev = pc_prev.points[p];

		mp_prev.dLdv_next = dt * mp.dLdx + mp.dLdv;
		mp_prev.dLdC_next = dt * mp.dLdF * (mp_prev.F + mp_prev.dFc).transpose() + mp.dLdC;
	}


	// Back prop G2P
	for (size_t p = 0; p < pc.points.size(); p++) {
		const MaterialPoint& mp_prev = pc_prev.points[p];


		auto nodes = grid.QueryPoint_CubicBSpline(mp_prev.x);

		for (size_t i = 0; i < nodes.size(); i++) {
			GridNode& node = nodes[i];

			Vec3 xg = node.x;
			Vec3 xp = mp_prev.x;
			Vec3 dgp = xg - xp;
			double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

			node.dLdv += mp_prev.dLdv_next * wgp + 3.0 / (dx * dx) * wgp * mp_prev.dLdC_next * dgp;
		}
	}

	// Back prop G_op
	for (size_t i = 0; i < grid.dim_x; i++) {
		for (size_t j = 0; j < grid.dim_y; j++) {
			for (size_t k = 0; k < grid.dim_z; k++) {
				GridNode& node = grid.nodes[i][j][k];

				if (node.m != 0.0) {
					node.dLdp = node.dLdv / node.m;

					node.dLdm = -1.0 / node.m * node.v.dot(node.dLdv);
				}
			}
		}
	}

	

	
	for (size_t p = 0; p < pc.points.size(); p++) {
		const MaterialPoint& mp = pc.points[p];
		MaterialPoint& mp_prev = pc_prev.points[p];

		// Back prop P2G
		auto nodes = grid.QueryPoint_CubicBSpline(mp_prev.x);

		for (size_t i = 0; i < nodes.size(); i++) {
			GridNode& node = nodes[i];

			Vec3 xg = node.x;
			Vec3 xp = mp_prev.x;
			Vec3 dgp = xg - xp;
			double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);
			Vec3 wgpGrad = -1.0 / dx * Vec3(
				CubicBSplineSlope(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx),
				CubicBSpline(dgp[0] / dx) * CubicBSplineSlope(dgp[1] / dx) * CubicBSpline(dgp[2] / dx),
				CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSplineSlope(dgp[2] / dx)
			);

			Mat3 G = -3.0 / (dx * dx) * dt * mp_prev.vol * mp_prev.P * (mp_prev.F + mp_prev.dFc).transpose() + mp_prev.m * mp_prev.C;

			// dL / dP
			mp_prev.dLdP -= wgp * 3.0 / (dx * dx) * dt * mp_prev.vol * node.dLdp * ((mp_prev.F + mp_prev.dFc).transpose() * dgp).transpose();

			// dL / dF
			// Note: this gradient also gets particle contributions, which will be added later
			mp_prev.dLdF -= wgp * 3.0 / (dx * dx) * dt * mp_prev.vol * dgp * (mp_prev.P.transpose() * node.dLdp).transpose();

			// dL / dC
			mp_prev.dLdC += wgp * mp_prev.m * node.dLdp * dgp.transpose();

			// dL / dm
			mp_prev.dLdm += wgp * node.dLdm;

			// dL / dx
			// Note: this gradient also gets particle contributions, which will be added later
			mp_prev.dLdx += mp_prev.m * node.dLdm * wgpGrad;
			mp_prev.dLdx += (wgpGrad * (mp_prev.m * mp_prev.v + G * dgp).transpose() - wgp * G.transpose()) * node.dLdp;
			mp_prev.dLdx += wgpGrad * node.v.transpose() * mp_prev.dLdv_next;
			Vec3 temp = InnerProduct(mp_prev.dLdC_next, node.v * dgp.transpose()) * wgpGrad;
			temp -= wgp * mp_prev.dLdC_next.transpose() * node.v;
			mp_prev.dLdx += 3.0 / (dx * dx) * temp;

			// dL / dv
			// mp_prev.dLdv += wgp * mp_prev.m * node.dLdp; // PREVIOUS GOOD ENOUGH VERSION
			mp_prev.dLdv += wgp * mp_prev.m * (1.0 - dt * drag) * node.dLdp; // CURRENT TRYING TO GET WORKING VERSION
		}


		// Back prop P_op_1
		mp_prev.dLdF += (Mat3::Identity() + dt * mp.C).transpose() * mp.dLdF;
		mp_prev.dLdF += d2_FCE_psi_dF2_mult_by_dF(mp_prev.F + mp_prev.dFc, mp_prev.lam, mp_prev.mu, mp_prev.dLdP);

		mp_prev.dLdx += mp.dLdx;
		mp_prev.dLdm += mp.dLdm;
	}
}
