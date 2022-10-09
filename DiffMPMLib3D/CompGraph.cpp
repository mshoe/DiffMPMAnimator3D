#include "pch.h"
#include "CompGraph.h"
#include "ForwardSimulation.h"
#include "BackPropagation.h"
#include "Interpolation.h"

void CompGraph::ComputeForwardPass(size_t start_layer)
{
	for (size_t i = start_layer; i < layers.size() - 1; i++)
	{
		ForwardTimeStep(
			*layers[i + (size_t)1].point_cloud,
			*layers[i].point_cloud,
			*layers[i].grid,
			dt, drag, f_ext);
	}
}

void CompGraph::ComputeBackwardPass(size_t control_layer)
{
	for (int i = (int)layers.size() - 2; i >= (int)control_layer; i--)
	{
		layers[i].grid->ResetGradients();
		layers[i].point_cloud->ResetGradients();

		Back_Timestep(layers[i + 1], layers[i], drag, dt);
	}
}

CompGraph::CompGraph(std::shared_ptr<PointCloud> initial_point_cloud, std::shared_ptr<Grid> grid, std::shared_ptr<const Grid> _target_grid)
{
	layers.clear();
	layers.resize(1);

	layers[0].point_cloud = std::make_shared<PointCloud>(*initial_point_cloud);
	layers[0].grid = std::make_shared<Grid>(*grid);

	target_grid = _target_grid;
}

void CompGraph::OptimizeDefGradControlSequence(
	// SIMULATION
	int num_steps, // number of timestepes, aka layers in the comp graph
	double _dt,
	double _drag,
	Vec3 _f_ext,

	// OPTIMIZATION
	int control_stride, // stride between control frames
	int max_gd_iters, // max gradient descent iterations per control frame
	int max_line_search_iters, // max line search iterations per control frame
	double initial_alpha, // initial gradient descent step size
	double gd_tol // tolerance factor used for determining when gradient descent has converged
)
{
	std::streamsize prev_precision = std::cout.precision(16);
	dt = _dt;
	drag = _drag;
	f_ext = _f_ext;


	// Initialize the computation graph
	layers.resize(num_steps);
	for (size_t i = 1; i < num_steps; i++) {
		layers[i].point_cloud = std::make_shared<PointCloud>(*layers.front().point_cloud);
		layers[i].grid = std::make_shared<Grid>(*layers.front().grid);
	}

	auto MassLoss = [this]
	{
		// Computes both the mass loss and the initial derivatives for backpropagation
		assert(target_grid->dim_x == grid->dim_x &&
			target_grid->dim_y == grid->dim_y &&
			target_grid->dim_z == grid->dim_z
		);

		PointCloud& point_cloud = *layers.back().point_cloud;
		Grid& grid = *layers.back().grid;
		double dx = grid.dx;
		P2G(point_cloud, grid, 0.0, 0.0);

		double loss = 0.0;
		for (size_t i = 0; i < target_grid->dim_x; i++) {
			for (size_t j = 0; j < target_grid->dim_y; j++) {
				for (size_t k = 0; k < target_grid->dim_z; k++) {
					// L
					double dm = grid.nodes[i][j][k].m - target_grid->nodes[i][j][k].m;
					loss += 0.5 * dm * dm;

					// dLdm
					grid.nodes[i][j][k].dLdm = dm;
				}
			}
		}

		// dLdx per particle
		for (size_t p = 0; p < point_cloud.points.size(); p++) {
			MaterialPoint& mp = point_cloud.points[p];
			mp.dLdx.setZero();

			auto nodes = grid.QueryPoint_CubicBSpline(mp.x);

			for (size_t i = 0; i < nodes.size(); i++) {
				const GridNode& node = nodes[i];

				Vec3 xg = node.x;
				Vec3 xp = mp.x;
				Vec3 dgp = xg - xp;
				double wgp = CubicBSpline(dgp[0] / dx) * CubicBSpline(dgp[1] / dx) * CubicBSpline(dgp[2] / dx);

				Vec3 wgpGrad = -1.0 / dx * Vec3(
					CubicBSplineSlope(dgp[0] / dx)	* CubicBSpline(dgp[1] / dx)			* CubicBSpline(dgp[2] / dx),
					CubicBSpline(dgp[0] / dx)		* CubicBSplineSlope(dgp[1] / dx)	* CubicBSpline(dgp[2] / dx),
					CubicBSpline(dgp[0] / dx)		* CubicBSpline(dgp[1] / dx)			* CubicBSplineSlope(dgp[2] / dx)
				);

				mp.dLdx += mp.m * node.dLdm * wgpGrad;
			}


			// other derivatives
			mp.dLdF.setZero();
			mp.dLdv.setZero();
			mp.dLdC.setZero();
		}

		return loss;
	};

	ComputeForwardPass(0);
	double initial_loss = MassLoss();
	std::cout << "initial loss = " << initial_loss << std::endl;
	ComputeBackwardPass(0);
	double initial_norm = layers.front().point_cloud->Compute_dLdF_Norm();
	std::cout << "initial norm = " << initial_norm << std::endl;
	double convergence_norm = gd_tol * initial_norm; // if norm is less than this, we say the timestep has converged
	std::cout << "convergence norm = " << convergence_norm << std::endl;

//#define FD 1
#ifdef FD
	/*******FINITE DIFFERENCES TEST********/

	MaterialPoint& mp = layers.front().point_cloud->points[0];

	std::cout << "analytic dLdF for particle 0:\n" << mp.dLdF << std::endl;


	Mat3 fd_dLdF;
	double delta = 1e-6;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			//MPMForwardSimulation(stcg, f_ext, drag, dt, control_timestep, false);
			double originalValue = mp.F(i, j);

			mp.F(i, j) = originalValue + delta;
			ComputeForwardPass(0);
			double l1 = MassLoss();

			mp.F(i, j) = originalValue - delta;
			ComputeForwardPass(0);
			double l2 = MassLoss();

			fd_dLdF(i, j) = (l1 - l2) / (2.0 * delta);
			mp.F(i, j) = originalValue;
		}
	}
	std::cout << "finite differences dLdF for particle 0:\n" << fd_dLdF << std::endl;

	double grad_diff = (fd_dLdF - mp.dLdF).norm();
	std::cout << "gradient difference = " << grad_diff << std::endl;

	double grad_diff_percent = 100.0 * grad_diff / fd_dLdF.norm();
	std::cout << "gradient difference percentage = " << grad_diff_percent << std::endl;
#else


	for (int control_timestep = 0; control_timestep < (int)layers.size() - 1; control_timestep += control_stride)
	{
		std::cout << "Optimizing for timestep: " << control_timestep << std::endl;

		double alpha = initial_alpha;
		for (int gd_iter = 0; gd_iter < max_gd_iters; gd_iter++)
		{
			std::cout << "gradient descent iteration = " << gd_iter << std::endl;

			// 1. compute forward pass from the control timestep
			ComputeForwardPass(control_timestep);

			// 2. get the current loss and final layer loss gradients
			double gd_loss = MassLoss();
			std::cout << "gd loss = " << gd_loss << std::endl;

			// 3. propagate loss gradients to control timestep
			ComputeBackwardPass(control_timestep);

			// 4. get the gradient norm
			double gradient_norm = layers.front().point_cloud->Compute_dLdF_Norm();

			// TODO: check for NaN dLdF

			// 5. Check if we have converged
			if (gradient_norm < gd_tol * initial_norm) {
				std::cout << "gradient norm = " << gradient_norm << " < " << convergence_norm << std::endl;
				std::cout << "control timestep " << control_timestep << " converged. Exiting gradient descent and moving to next control timestep." << std::endl;
				break;
			}

			// 6. Line search to determine a step size where the loss function decreases
			double ls_loss = gd_loss;
			bool ls_loss_decrease_found = false;
			for (int ls_iter = 0; ls_iter < max_line_search_iters; ls_iter++) 
			{

				layers[control_timestep].point_cloud->Descend_dLdF(alpha, gradient_norm);

				ComputeForwardPass(control_timestep);
				ls_loss = MassLoss();

				if (ls_loss < gd_loss) 
				{
					ls_loss_decrease_found = true;
					std::cout << "line search decrease found at ls_iter = " << ls_iter << ", loss = " << ls_loss << std::endl;
					break;
				}

				layers[control_timestep].point_cloud->Descend_dLdF(-alpha, gradient_norm);

				alpha /= 2.0;
			}


			if (!ls_loss_decrease_found) {
				std::cout << "Line search unable to find a decrease in timestep (" << control_timestep << 
					"), gd iteration (" << gd_iter << ")" << std::endl;
				std::cout << "Exiting gradient descent and moving to next control timestep." << std::endl;
				break;
			}
		}
	}

	double final_loss = MassLoss();
	std::cout << "Initial loss = " << initial_loss << std::endl;
	std::cout << "Final loss = " << final_loss << std::endl;

	if (final_loss >= initial_loss) {
		std::cout << "Unable to find a loss decrease." << std::endl;
	}

	std::cout << "Deformation gradient control sequence optimization finished." << std::endl;
#endif
	std::cout.precision(prev_precision);
}