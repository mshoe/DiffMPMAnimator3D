#include "pch.h"
#include "SphereUnionSurfacing.h"
#include "Grid.h"
#include "igl/marching_cubes.h"
#include <chrono>

void DiffMPMLib3D::SphereUnionMarchingCubesSurfaceFromPointCloud(const std::vector<Vec3>& _points, double radius, double grid_dx, double iso_mass, int blur_iterations, Vec3 grid_min_point, Vec3 grid_max_point, Eigen::MatrixXd& mcV, Eigen::MatrixXi& mcF)
{
    std::cout << "generating marching cubes surface..." << std::endl;
    auto begin_clock = std::chrono::steady_clock::now();

    Eigen::Vector3i grid_dims;
    grid_dims.resize(3);
    for (int i = 0; i < 3; i++) {
        grid_dims[i] = (int)std::ceil((grid_max_point[0] - grid_min_point[0]) / grid_dx);
    }
    std::vector<std::vector<std::vector<double>>> grid;
    grid.resize(grid_dims[0]);

    for (int i = 0; i < grid_dims[0]; i++) {
        grid[i].resize(grid_dims[1]);
        for (int j = 0; j < grid_dims[1]; j++) {
            grid[i][j].resize(grid_dims[2]);
            for (int k = 0; k < grid_dims[2]; k++) {
                grid[i][j][k] = 0.0;
            }
        }
    }

    Vec3 radius_vec = Vec3(radius, radius, radius);

    // for each point, find the nodes that intersect with its "sphere of influence", then set the mass at those nodes to 1
    for (size_t p = 0; p < _points.size(); p++) {
        Vec3 point = _points[p];
        Vec3 aabb_min_point = point - radius_vec;
        Vec3 aabb_max_point = point + radius_vec;

        // Get the indices that overlap with this aabb
        Eigen::Vector3i min_index;
        for (int i = 0; i < 3; i++) {
            min_index[i] = (int)std::ceil((aabb_min_point[i] - grid_min_point[i]) / grid_dx);
        }

        Eigen::Vector3i max_index;
        for (int i = 0; i < 3; i++) {
            max_index[i] = (int)std::ceil((aabb_max_point[i] - grid_min_point[i]) / grid_dx);
        }

        // clip min and max indices to be in bounds
        for (int i = 0; i < 3; i++) {
            min_index[i] = std::max(0, min_index[i]);
            max_index[i] = std::min(grid_dims[i] - 1, max_index[i]);
        }
        /*std::cout << "min_index = " << min_index.transpose() << std::endl;
        std::cout << "max_index = " << min_index.transpose() << std::endl;*/
        // now check if the node at each index is inside the sphere.
        // If yes, then set mass to 1
        for (int i = min_index[0]; i <= max_index[0]; i++) {
            for (int j = min_index[1]; j <= max_index[1]; j++) {
                for (int k = min_index[2]; k <= max_index[2]; k++) {
                    Vec3 curr_index = Vec3(i, j, k);
                    Vec3 node_pos = grid_min_point + grid_dx * curr_index;
                    /*std::cout << "node_pos = " << node_pos.transpose() << std::endl;
                    std::cout << "point = " << point.transpose() << std::endl;*/

                    if ((node_pos - point).squaredNorm() <= radius * radius) {
                        
                        // distance to sphere center
                        grid[i][j][k] = std::max(grid[i][j][k], radius - (node_pos - point).norm());
                    }
                }
            }
        }
    }

    // Now do some blurring/smoothing
    auto blur_grid = grid;
    for (int blur_iteration = 0; blur_iteration < blur_iterations; blur_iteration++) {
        for (int i = 1; i < grid_dims[0] - 1; i++) {
            for (int j = 1; j < grid_dims[1] - 1; j++) {
                for (int k = 1; k < grid_dims[2] - 1; k++) {

                    double avg_value = 0.0;
                    for (int ii = -1; ii <= 1; ii++) {
                        for (int jj = -1; jj <= 1; jj++) {
                            for (int kk = -1; kk <= 1; kk++) {
                                avg_value += grid[i+ii][j+jj][k+kk];
                            }
                        }
                    }
                    avg_value /= 9.0;

                    blur_grid[i][j][k] = avg_value;
                }
            }
        }

        grid = blur_grid;
    }

    // Now construct GV and Gf
    Eigen::MatrixXd GV;
    Eigen::VectorXd Gf;
    GV.resize(grid_dims[0] * grid_dims[1] * grid_dims[2], 3);
    Gf.resize(grid_dims[0] * grid_dims[1] * grid_dims[2]);
    size_t ind = 0;
    for (size_t i = 0; i < grid_dims[0]; i++) {
        for (size_t j = 0; j < grid_dims[1]; j++) {
            for (size_t k = 0; k < grid_dims[2]; k++) {
                Vec3 curr_index = Vec3(i, j, k);
                Vec3 node_pos = grid_min_point + grid_dx * curr_index;
                GV.row(ind) = node_pos;

                Gf(ind) = grid[i][j][k];
                /*if (Gf(ind) != 0.0) {
                    std::cout << curr_index.transpose() << std::endl;
                }*/
                ind++;
            }
        }
    }
    
    // Now construct marching cubes meshes
    igl::marching_cubes(Gf, GV, grid_dims[0], grid_dims[1], grid_dims[2], iso_mass, mcV, mcF);

    auto end_clock = std::chrono::steady_clock::now();
    std::cout << "Finished generating marching cubes surface in " << std::chrono::duration_cast<std::chrono::seconds>(end_clock - begin_clock).count() << std::endl;
}
