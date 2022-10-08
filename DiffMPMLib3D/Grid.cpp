#include "pch.h"
#include "Grid.h"

Grid::Grid(int _dim_x, int _dim_y, int _dim_z, double _dx, Vec3 _min_point)
	: dim_x(_dim_x), dim_y(_dim_y), dim_z(_dim_z), dx(_dx), min_point(_min_point)
{
	nodes = std::vector<std::vector<std::vector<GridNode>>>();
	nodes.resize(dim_x);

	for (int i = 0; i < dim_x; i++) {
		nodes[i].resize(dim_y);
		for (int j = 0; j < dim_y; j++) {
			nodes[i][j].resize(dim_z);
			for (int k = 0; k < dim_z; k++) {
				// potentially not worth storing this x...
				nodes[i][j][k].x = min_point + Vec3(i, j, k) * dx;
			}
		}
	}
}

std::vector<std::reference_wrapper<GridNode>> Grid::QueryPoint_CubicBSpline(Vec3 point)
{
	/* 
	* Returns all the nodes which are within interpolation range of the point position
	*/

	std::vector<std::reference_wrapper<GridNode>> ret;

	Vec3 relative_point = point - min_point;

	int bot_left_index[3];
	for (size_t i = 0; i < 3; i++) {
		bot_left_index[i] = (int)std::floor(relative_point[i] / dx) - 1;
	}

	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++) {
			for (int k = 0; k <= 3; k++) {
				int index[3] = {
					bot_left_index[0] + i,
					bot_left_index[1] + j,
					bot_left_index[2] + k
				};

				// check if this is a valid node
				if (0 <= index[0] && index[0] < dim_x &&
					0 <= index[1] && index[1] < dim_y &&
					0 <= index[2] && index[2] < dim_z)
				{
					ret.push_back(nodes[index[0]][index[1]][index[2]]);
				}
			}
		}
	}

	return ret;
}

std::vector<Vec3> Grid::GetNodePositions() const
{
	std::vector<Vec3> ret;
	for (int i = 0; i < dim_x; i++) {
		for (int j = 0; j < dim_y; j++) {
			for (int k = 0; k < dim_z; k++) {
				ret.push_back(nodes[i][j][k].x);
			}
		}
	}
	return ret;
}