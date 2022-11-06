#include "pch.h"
#include "MultiThreadForwardSimulation.h"
#include "ForwardSimulation.h"
#include <thread>
#include <functional>

void DiffMPMLib3D::MultiThreadMPM::ForwardTimeStep(PointCloud& next_point_cloud, PointCloud& curr_point_cloud, Grid& main_grid,
	std::vector<std::shared_ptr<Grid>>& proxy_grids, double dt, double drag, Vec3 f_ext)
{
	size_t num_threads = proxy_grids.size();
	/*bool good_num_threads = (num_threads == 1 || num_threads == 2 || num_threads == 4 || num_threads == 8);
	if (!good_num_threads)
	{
		std::cout << "error: number of threads not a power of 2" << std::endl;
		return;
	}*/

	/**** PARTICLE THREAD SPLITTING ****/
	// split the particles into num_threads sections (This could probably be done outside, but probably isn't much overhead here anyways
	std::vector<std::pair<size_t, size_t>> p_sections; // first is start index, second is number of particles
	
	size_t num_points = curr_point_cloud.points.size();
	size_t min_num_points_per_thread = num_points / num_threads;
	size_t num_remainder_points = num_points % num_threads; // just add this to the first thread

	size_t next_start_ind = min_num_points_per_thread + num_remainder_points;
	p_sections.push_back({ 0, next_start_ind });

	for (size_t i = 1; i < num_threads; i++) {
		p_sections.push_back({ next_start_ind, min_num_points_per_thread });
		next_start_ind += min_num_points_per_thread;
	}


	/**** GRID THREAD SPLITTING ****/
	struct GridThreadInfo
	{
		size_t grid_block_start_i = 0;
		size_t grid_block_start_j = 0;
		size_t grid_block_start_k = 0;

		size_t grid_block_size_i = 0;
		size_t grid_block_size_j = 0;
		size_t grid_block_size_k = 0;
	};
	std::vector<GridThreadInfo> grid_section;
	//size_t min_nodes_i_per_thread = (size_t)main_grid.dim_x / num_threads;
	//size_t min_nodes_j_per_thread = (size_t)main_grid.dim_y / num_threads;
	size_t min_nodes_k_per_thread = (size_t)main_grid.dim_z / num_threads;
	//size_t num_remainder_nodes_i = (size_t)main_grid.dim_x % num_threads;
	//size_t num_remainder_nodes_j = (size_t)main_grid.dim_y % num_threads;
	size_t num_remainder_nodes_k = (size_t)main_grid.dim_z % num_threads;

	//size_t next_start_ind_i = min_nodes_i_per_thread + num_remainder_nodes_i;
	//size_t next_start_ind_j = min_nodes_j_per_thread + num_remainder_nodes_j;
	size_t next_start_ind_k = min_nodes_k_per_thread + num_remainder_nodes_k;
	grid_section.push_back({ 0, 0, 0, (size_t)main_grid.dim_x, (size_t)main_grid.dim_y, next_start_ind_k});
	for (size_t i = 1; i < num_threads; i++) {
		grid_section.push_back({ 0, 0, next_start_ind_k, (size_t)main_grid.dim_x, (size_t)main_grid.dim_z, min_nodes_k_per_thread});
		//next_start_ind_i += min_nodes_i_per_thread;
		//next_start_ind_j += min_nodes_j_per_thread;
		next_start_ind_k += min_nodes_k_per_thread;
	}

	//std::cout << "num threads = " << num_threads << std::endl;
	//std::cout << "min nodes i per thread = " << min_nodes_i_per_thread << std::endl;
	//std::cout << "num remainder nodes i = " << num_remainder_nodes_i << std::endl;
	/*std::cout << "grid dim x = " << main_grid.dim_x << std::endl;
	std::cout << "f_ext = " << f_ext.transpose() << std::endl;
	std::cout << "dt = " << dt << std::endl;
	std::cout << "min num points per thread = " << min_num_points_per_thread << std::endl;
	std::cout << "num points = " << num_points << std::endl;
	std::cout << "num points (calculated) = " << p_sections.back().first + p_sections.back().second << std::endl;*/

	std::vector<std::thread> parallel_threads;

	/** P OP 1 **/
	for (size_t i = 0; i < num_threads; i++) 
	{
		parallel_threads.push_back(std::thread(MultiThreadMPM::Multi_P_op_1, std::ref(curr_point_cloud.points), p_sections[i].first, p_sections[i].second));
	}
	for (size_t i = 0; i < num_threads; i++) 
	{
		parallel_threads[i].join();
	}
	parallel_threads.clear();


	/** RESET PROXY GRIDS **/
	for (size_t i = 0; i < num_threads; i++) 
	{
		parallel_threads.push_back(std::thread(SingleThreadMPM::G_Reset, std::ref(*proxy_grids[i])));
	}
	for (size_t i = 0; i < num_threads; i++)
	{
		parallel_threads[i].join();
	}
	parallel_threads.clear();


	/** P 2 G **/
	for (size_t i = 0; i < num_threads; i++) 
	{
		parallel_threads.push_back(std::thread(MultiThreadMPM::Multi_P2G,
			std::cref(curr_point_cloud.points), std::ref(*proxy_grids[i]),
			p_sections[i].first, p_sections[i].second,
			dt, drag));
	}
	for (size_t i = 0; i < num_threads; i++) {
		parallel_threads[i].join();
	}
	parallel_threads.clear();

	/** Grid Operation **/
	for (size_t i = 0; i < num_threads; i++)
	{
		parallel_threads.push_back(std::thread(MultiThreadMPM::Multi_G_op, std::cref(proxy_grids), std::ref(main_grid),
			grid_section[i].grid_block_start_i, grid_section[i].grid_block_start_j, grid_section[i].grid_block_start_k,
			grid_section[i].grid_block_size_i, grid_section[i].grid_block_size_j, grid_section[i].grid_block_size_k,
			dt, f_ext));
	}
	for (size_t i = 0; i < num_threads; i++)
	{
		parallel_threads[i].join();
	}
	parallel_threads.clear();

	/** G 2 P **/
	for (size_t i = 0; i < num_threads; i++)
	{
		parallel_threads.push_back(std::thread(MultiThreadMPM::Multi_G2P,
			std::ref(next_point_cloud.points), std::cref(curr_point_cloud.points), std::cref(main_grid),
			p_sections[i].first, p_sections[i].second));
	}
	for (size_t i = 0; i < num_threads; i++) {
		parallel_threads[i].join();
	}
	parallel_threads.clear();

	/** P OP 2 **/
	for (size_t i = 0; i < num_threads; i++)
	{
		parallel_threads.push_back(std::thread(MultiThreadMPM::Multi_P_op_2,
			std::ref(next_point_cloud.points), std::cref(curr_point_cloud.points),
			p_sections[i].first, p_sections[i].second,
			dt));
	}
	for (size_t i = 0; i < num_threads; i++) {
		parallel_threads[i].join();
	}
	parallel_threads.clear();
}

void DiffMPMLib3D::MultiThreadMPM::Multi_P_op_1(std::vector<MaterialPoint>& points, size_t p_ind_start, size_t num_p)
{
	for (size_t p = p_ind_start; p < p_ind_start + num_p; p++) {

		// MP_op_1
		MaterialPoint& mp = points[p];

		SingleThreadMPM::SingleParticle_op_1(mp);
	}
}


void DiffMPMLib3D::MultiThreadMPM::Multi_P2G(const std::vector<MaterialPoint>& points, Grid& proxy_grid, size_t p_ind_start, size_t num_p, double dt, double drag)
{
	for (size_t p = p_ind_start; p < p_ind_start + num_p; p++) {

		// MP_op_1
		const MaterialPoint& mp = points[p];

		SingleThreadMPM::SingleParticle_to_grid(mp, proxy_grid, dt, drag);
	}
}

void DiffMPMLib3D::MultiThreadMPM::Multi_G_op(const std::vector<std::shared_ptr<Grid>>& proxy_grids, Grid& main_grid,
	size_t block_start_i, size_t block_start_j, size_t block_start_k,
	size_t block_size_i, size_t block_size_j, size_t block_size_k,
	double dt, Vec3 f_ext)
{
	for (size_t i = block_start_i; i < block_start_i + block_size_i; i++) {
		for (size_t j = block_start_j; j < block_start_j + block_size_j; j++) {
			for (size_t k = block_start_k; k < block_start_k + block_size_k; k++) {

				// RESET THE MAIN GRID
				main_grid.nodes[i][j][k].m = 0.0;
				main_grid.nodes[i][j][k].v.setZero();
				main_grid.nodes[i][j][k].p.setZero();

				// REDUCE THE PROXY GRIDS TO THE MAIN GRID
				for (size_t l = 0; l < proxy_grids.size(); l++) {
					main_grid.nodes[i][j][k].m += proxy_grids[l]->nodes[i][j][k].m;
					main_grid.nodes[i][j][k].p += proxy_grids[l]->nodes[i][j][k].p;
				}
				SingleThreadMPM::SingleNode_op(main_grid.nodes[i][j][k], dt, f_ext);
			}
		}
	}
}

void DiffMPMLib3D::MultiThreadMPM::Multi_G2P(std::vector<MaterialPoint>& next_points, const std::vector<MaterialPoint>& curr_points, const Grid& main_grid,
	size_t p_ind_start, size_t num_p)
{
	for (size_t p = p_ind_start; p < p_ind_start + num_p; p++) 
	{
		const MaterialPoint& curr_mp = curr_points[p];
		MaterialPoint& next_mp = next_points[p];

		SingleThreadMPM::Grid_to_SingleParticle(next_mp, curr_mp, main_grid);
	}
}

void DiffMPMLib3D::MultiThreadMPM::Multi_P_op_2(std::vector<MaterialPoint>& next_points, const std::vector<MaterialPoint>& curr_points, size_t p_ind_start, size_t num_p, double dt)
{
	for (size_t p = p_ind_start; p < p_ind_start + num_p; p++)
	{
		const MaterialPoint& curr_mp = curr_points[p];
		MaterialPoint& next_mp = next_points[p];

		SingleThreadMPM::SingleParticle_op_2(next_mp, curr_mp, dt);
	}
}
