#pragma once

#include "PointCloud.h"
#include "Grid.h"
#include "CompGraph.h"
#include "polyscope/point_cloud.h"



struct SceneInput
{
	std::string mpm_input_mesh_path = "experiments\\004 sphere to bunny\\input\\icosphere.obj";
	std::string mpm_target_mesh_path = "experiments\\004 sphere to bunny\\input\\bunny.obj";


	double grid_dx = 1.0; // TODO: make this 1.0 and use scaled up meshes
	int points_per_cell_cuberoot = 2;

	/*Vec3 grid_min_point = Vec3(-8.0, -8.0, -8.0);
	Vec3 grid_max_point = Vec3(8.0, 8.0, 8.0);*/
	Vec3 grid_min_point = Vec3(-8.0, -8.0, -8.0) * 2.0;
	Vec3 grid_max_point = Vec3(8.0, 8.0, 8.0) * 2.0;

	double lam = 58333.0;
	double mu = 38888.9;
	//double lam = 2500.0;
	//double mu = 100.0;
	double p_density = 40.0;// *5.0;

	double dt = 1.0 / 120.0;// / 100.0;
	double drag = 0.5;
	//Vec3 f_ext = Vec3(0.0, 0.0, -9.81);
	Vec3 f_ext = Vec3(0.0, 0.0, 0.0);
};

bool LoadMPMPointCloudFromObj(
	std::string obj_path,
	std::shared_ptr<PointCloud>& mpm_point_cloud,
	double point_dx,
	double density,
	double lam,
	double mu
);

bool LoadScene(
	const SceneInput& scene_input,
	std::shared_ptr<PointCloud>& mpm_point_cloud, 
	std::shared_ptr<Grid>& mpm_grid,
	polyscope::PointCloud** polyscope_point_cloud,
	polyscope::PointCloud** polyscope_grid
);

bool LoadCompGraph(
	const SceneInput& scene_input,
	std::shared_ptr<CompGraph>& comp_graph,
	polyscope::PointCloud** polyscope_input_point_cloud,
	polyscope::PointCloud** polyscope_target_point_cloud,
	polyscope::PointCloud** polyscope_grid
);