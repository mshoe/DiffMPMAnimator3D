// PointCloudMesher.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include <iostream>
#include <fstream>

#include "ForwardSimulation.h"
using namespace DiffMPMLib3D;
//#include "igl/dual_contouring.h"
#include "Eigen/Dense"

//#include "libigl_tutorial_715.h"

#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <igl/parallel_for.h>
#include <igl/writeOBJ.h>
#include <igl/per_vertex_normals.h>

int num_objs = 600;

std::string input_path = "experiments\\003 sphere to spot\\output\\";
std::string output_path = "experiments\\003 sphere to spot\\meshes\\";

polyscope::PointCloud* ps_point_cloud = nullptr;
polyscope::SurfaceMesh* ps_surface_mesh = nullptr;
std::vector<Eigen::Vector3d> points;

unsigned GetNumberOfDigits(unsigned i)
{
    return i > 0 ? (int)log10((double)i) + 1 : 1;
}

std::string LeadingZerosNumberStr(int number, int num_digits_of_string)
{
    std::string ret = "";
    int num_digits_of_number = GetNumberOfDigits(number);
    for (unsigned int j = 0; j < num_digits_of_string - num_digits_of_number; j++) {
        ret += "0";
    }
    ret += std::to_string(number);
    return ret;
}


void LoadPointCloud(std::string path_to_points, std::vector<Vec3>& _points)
{
    _points.clear();
    std::ifstream ifs;
    ifs.open(path_to_points);


    if (ifs.good()) 
    {
        std::string junk;
        double x, y, z;

        while (ifs >> junk >> x >> y >> z)
        {
            Eigen::Vector3d _points(x, y, z);
            points.push_back(_points);
        }
    }
    else {
        std::cout << "couldn't open file" << std::endl;
    }
}

void MarchingCubesPointCloud(const std::vector<Vec3>& _points, double iso_mass, double grid_dx, Vec3 grid_min_point, Vec3 grid_max_point,
    Eigen::MatrixXd& mcV, Eigen::MatrixXi& mcF)
{
    std::cout << "generating marching cubes surface..." << std::endl;
    auto begin_clock = std::chrono::steady_clock::now();
    int grid_dims[3];
    for (int i = 0; i < 3; i++) {
        grid_dims[i] = std::ceil((grid_max_point[0] - grid_min_point[0]) / grid_dx);
    }
    auto mpm_grid = std::make_shared<Grid>(grid_dims[0], grid_dims[1], grid_dims[2], grid_dx, grid_min_point);

    SingleThreadMPM::P2G_Mass(points, *mpm_grid, 1.0);

    double quality = 4.0;
    auto high_res_mpm_grid = std::make_shared<Grid>(grid_dims[0] * quality, grid_dims[1] * quality, grid_dims[2] * quality, grid_dx / quality, grid_min_point);

    SingleThreadMPM::G2G_Mass(*mpm_grid, *high_res_mpm_grid);

    

    Eigen::MatrixXd GV; // location of each grid node
    Eigen::VectorXd Gf;
    high_res_mpm_grid->GetMassSDF(GV, Gf);

    igl::marching_cubes(Gf, GV, high_res_mpm_grid->dim_x, high_res_mpm_grid->dim_y, high_res_mpm_grid->dim_z, iso_mass, mcV, mcF);

    auto end_clock = std::chrono::steady_clock::now();
    std::cout << "Finished generating marching cubes surface in " << std::chrono::duration_cast<std::chrono::seconds>(end_clock - begin_clock).count() << std::endl;
}

void mesherCallback()
{
    ImGui::PushItemWidth(100);

    static double grid_dx = 1.0;
    static Vec3 grid_min_point = Vec3(-16.0, -16.0, -16.0);
    static Vec3 grid_max_point = Vec3(16.0, 16.0, 16.0);

    static double iso_mass = 0.5;
    ImGui::InputDouble("grid dx", &grid_dx);
    ImGui::InputDouble("iso mass", &iso_mass);

    if (ImGui::Button("Generate meshes from point clouds"))
    {
        polyscope::SurfaceMesh* temp_surf = nullptr;
        for (int i = 0; i < num_objs; i++)
        {
            
            std::string pc_path = input_path + "mpm_points_" + LeadingZerosNumberStr(i, 6) + ".obj";
            std::string mesh_path = output_path + "mpm_isosurface_" + LeadingZerosNumberStr(i, 6) + ".obj";
            std::string ss_path = output_path + "screenshot_" + LeadingZerosNumberStr(i, 6) + ".png";

            std::cout << "loading " << pc_path << std::endl;

            LoadPointCloud(pc_path, points);
            Eigen::MatrixXd mcV;
            Eigen::MatrixXi mcF;
            MarchingCubesPointCloud(points, iso_mass, grid_dx, grid_min_point, grid_max_point, mcV, mcF);
            


            //std::cout << "writing to " << mesh_path << std::endl;
            //igl::writeOBJ(mesh_path, mcV, mcF);

            if (temp_surf)
                polyscope::removeSurfaceMesh("temp_surf");
            temp_surf = polyscope::registerSurfaceMesh("temp_surf", mcV, mcF);
            temp_surf->setSurfaceColor(glm::vec3(1.0, 0.5, 0.5));
            temp_surf->setSmoothShade(true);

            polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

            polyscope::screenshot(ss_path, false);
        }
    }

    if (ImGui::Button("Generate screenshots from point clouds"))
    {
        for (int i = 0; i < num_objs; i++)
        {

            std::string pc_path = input_path + "mpm_points_" + LeadingZerosNumberStr(i, 6) + ".obj";
            std::string ss_path = output_path + "screenshot_pointcloud_" + LeadingZerosNumberStr(i, 6) + ".png";

            std::cout << "loading " << pc_path << std::endl;

            LoadPointCloud(pc_path, points);
            auto pc = polyscope::registerPointCloud("curr cloud", points);
            
            pc->setPointRadius(0.02683);
            pc->setPointColor(glm::vec3(0.110, 0.388, 0.890));

            polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

            polyscope::screenshot(ss_path, false);

            polyscope::removePointCloud("curr cloud");
        }
    }


    if (ImGui::Button("load point cloud .obj"))
    {
        LoadPointCloud(input_path + "mpm_points_000599.obj", points);

        if (ps_point_cloud != nullptr) {
            polyscope::removePointCloud("point cloud");
        }
        ps_point_cloud = polyscope::registerPointCloud("point cloud", points);
    }


    if (!points.empty()) 
    {
        
        
        //ImGui::InputFloat3()
        if (ImGui::Button("marching cubes point cloud"))
        {
            if (ps_surface_mesh != nullptr)
                polyscope::removeStructure("surface mesh");

            Eigen::MatrixXd mcV;
            Eigen::MatrixXi mcF;
            MarchingCubesPointCloud(points, iso_mass, grid_dx, grid_min_point, grid_max_point, mcV, mcF);
            Eigen::MatrixXd N;
            igl::per_vertex_normals(mcV, mcF, N);

            ps_surface_mesh = polyscope::registerSurfaceMesh("surface mesh", mcV, mcF);
            //ps_surface_mesh->addVect
        }
    }

    ImGui::PopItemWidth();
}

int main()
{
    // Initialize polyscope, creating graphics contexts and constructing a window.
    // Should be called exactly once.
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::init();

    /*
    * build visualizations, here or in distant code
    *
    */


    // Specify the callback
    polyscope::state::userCallback = mesherCallback;

    // Pass control flow to polyscope, displaying the interactive window.
    // Function will return when user closes the window.
    polyscope::show();

    //libigl_tutorial_715();
}

