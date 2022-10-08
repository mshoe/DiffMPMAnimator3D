#pragma once
#include "pch.h"
#include "PointCloud.h"
#include "Grid.h"

struct CompGraphLayer
{
	std::shared_ptr<PointCloud> point_cloud = nullptr;
	std::shared_ptr<Grid> grid = nullptr;
};

class CompGraph
{
public:
	CompGraph(std::shared_ptr<PointCloud> initial_point_cloud, std::shared_ptr<Grid> grid, std::shared_ptr<const Grid> _target_grid);


	void OptimizeDefGradControlSequence(
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
	);


	void ComputeForwardPass(size_t start_layer);
	void ComputeBackwardPass();

	std::vector<CompGraphLayer> layers;
	
private:
	std::shared_ptr<const Grid> target_grid;
	Vec3 f_ext = Vec3::Zero();
	double dt = 1.0 / 120.0;
	double drag = 0.5;
};