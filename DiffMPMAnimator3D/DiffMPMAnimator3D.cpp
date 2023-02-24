// DiffMPMAnimator3D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"



// C++ includes
#include <iostream>
#include <sstream>
#include <chrono>

// DiffMPMLib3D
#include "MPMPolyscope.h"
#include "ForwardSimulation.h"
#include "MultiThreadForwardSimulation.h"
#include "Elasticity.h"
using namespace DiffMPMLib3D;

#include "cereal/archives/json.hpp"
#include "MXImGuiTools.h"

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

// GLOBAL VARS
OptInput opt_input;

// FOR REAL TIME
std::shared_ptr<PointCloud> mpm_point_cloud = nullptr;
std::shared_ptr<Grid> mpm_grid = nullptr;
polyscope::PointCloud* ps_point_cloud = nullptr;
polyscope::PointCloud* ps_grid = nullptr;

// FOR OPTIMIZATION
polyscope::PointCloud* ps_target_point_cloud = nullptr;
std::shared_ptr<CompGraph> comp_graph = nullptr;

// FOR POST-PROCESSING VISUALIZATION


void RealTimeMPMImGUI() 
{
    double screenshot_interval = 1.0 / 30.0;

    static double time = 0.0;
    static int num_threads = 4;
    static std::vector<std::shared_ptr<Grid>> proxy_grids;
    if (ImGui::InputInt("number of threads", &num_threads))
    {
        if (num_threads < 1) num_threads = 1;
        if (num_threads > 8) num_threads = 8;
    }
    static bool scene_loaded = false;
    if (scene_loaded)
        ImGui::BeginDisabled();
    if (ImGui::Button("Load Scene"))
    {
        LoadScene(opt_input, mpm_point_cloud, mpm_grid, &ps_point_cloud, &ps_grid);
        proxy_grids.clear();
        proxy_grids.resize(num_threads);
        for (int i = 0; i < num_threads; i++) {
            proxy_grids[i] = std::make_shared<Grid>(*mpm_grid);
        }
        scene_loaded = true;
    }
    if (scene_loaded)
        ImGui::EndDisabled();

    if (!scene_loaded)
    {
        ImGui::PopItemWidth();
        return;
    }

    static double& gravity = opt_input.f_ext[2];
    ImGui::InputDouble("gravity", &gravity);

    static bool multi_threaded = false;
    ImGui::Checkbox("multithreaded", &multi_threaded);
    if (ImGui::Button("timestep")) {
        if (!multi_threaded) {
            SingleThreadMPM::ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, opt_input.dt, opt_input.drag, opt_input.f_ext);
        }
        else {
            MultiThreadMPM::ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, proxy_grids, opt_input.dt, opt_input.drag, opt_input.f_ext);
        }
        time += opt_input.dt;

        auto point_positions = mpm_point_cloud->GetPointPositions();
        ps_point_cloud->updatePointPositions(point_positions);

        auto node_masses = mpm_grid->GetNodeMasses();
        ps_grid->addScalarQuantity("masses", node_masses);

        auto node_velocities = mpm_grid->GetNodeVelocities();
        ps_grid->addVectorQuantity("velocities", node_velocities);
    }

    static int num_timesteps = 100;
    ImGui::InputInt("Num timesteps", &num_timesteps);
    if (ImGui::Button("X timesteps"))
    {
        auto begin_clock = std::chrono::steady_clock::now();
        for (size_t i = 0; i < num_timesteps; i++) {
            if (!multi_threaded) {
                SingleThreadMPM::ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, opt_input.dt, opt_input.drag, opt_input.f_ext);
            }
            else {
                MultiThreadMPM::ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, proxy_grids, opt_input.dt, opt_input.drag, opt_input.f_ext);
            }
            time += opt_input.dt;
        }
        auto end_clock = std::chrono::steady_clock::now();
        std::cout << "Compute took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_clock - begin_clock).count() << " milliseconds." << std::endl;

        auto point_positions = mpm_point_cloud->GetPointPositions();
        ps_point_cloud->updatePointPositions(point_positions);

        auto node_masses = mpm_grid->GetNodeMasses();
        ps_grid->addScalarQuantity("masses", node_masses);

        auto node_velocities = mpm_grid->GetNodeVelocities();
        ps_grid->addVectorQuantity("velocities", node_velocities);
    }
    

    if (ImGui::Button("screenshot")) {
        polyscope::screenshot();
    }

    if (ImGui::Button("export point cloud points"))
    {
        mpm_point_cloud->WriteToOBJ("points.obj");
    }

    
    static bool playing = false;
    ImGui::Checkbox("playing", &playing);
    if (playing) {
        if (!multi_threaded) {
            SingleThreadMPM::ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, opt_input.dt, opt_input.drag, opt_input.f_ext);
        }
        else {
            MultiThreadMPM::ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, proxy_grids, opt_input.dt, opt_input.drag, opt_input.f_ext);
        }
        time += opt_input.dt;

        
        auto point_positions = mpm_point_cloud->GetPointPositions();
        ps_point_cloud->updatePointPositions(point_positions);
    }


    static int point_index = 0;
    ImGui::InputInt("point index", &point_index);

    if (point_index < 0) point_index = 0;
    if (point_index >= mpm_point_cloud->points.size()) point_index = mpm_point_cloud->points.size() - 1;

    const MaterialPoint& mp = mpm_point_cloud->points[point_index];

    std::stringstream ss;

    ss.str(std::string());
    ss << "time = " << time;
    ImGui::Text(ss.str().c_str());


    if (ImGui::TreeNode("mp properties")) {
        ss.str(std::string());
        ss << mp.x.transpose();
        ImGui::Text(std::string("x = " + ss.str()).c_str());

        ss.str(std::string());
        ss << mp.v.transpose();
        ImGui::Text(std::string("v = " + ss.str()).c_str());

        ss.str(std::string());
        ss << mp.F;
        ImGui::Text(std::string("F = \n" + ss.str()).c_str());

        ss.str(std::string());
        ss << mp.P;
        ImGui::Text(std::string("P = \n" + ss.str()).c_str());

        ss.str(std::string());
        ss << mp.C;
        ImGui::Text(std::string("C = \n" + ss.str()).c_str());



        ss.str(std::string());
        ss << mp.m;
        ImGui::Text(std::string("m = " + ss.str()).c_str());

        ss.str(std::string());
        ss << mp.vol;
        ImGui::Text(std::string("vol = " + ss.str()).c_str());

        ImGui::TreePop();
    }
}


void CheckGradients()
{
    static double lam = 58333.0;
    static double mu = 38888.9;
    if (ImGui::Button("Print gradients"))
    {
        // ROTATION MATRIX
        double a = 45.0 * 0.0174533;
        Mat3 F;
        F << 1, 0, 0, 0, cos(a), -sin(a), 0, sin(a), cos(a);

        //Mat3 F;
        //F.setIdentity();
        F *= 0.8;

        Tensor3x3x3x3 dJFit_dF = d_JFit_dF_FD(F);
        std::cout << "d(J * inv(F.T)/dF finite differences:\n";
        std::cout << ToString(dJFit_dF) << std::endl << std::endl;


        Tensor3x3x3x3 dP_dF_FD = d2_FCE_psi_dF2_FD(F, lam, mu);
        std::cout << "dP/dF finite differences:\n";
        std::cout << ToString(dP_dF_FD) << std::endl << std::endl;


        Tensor3x3x3x3 dP_dF_mult_trick = d2_FCE_psi_dF2_mult_trick(F, lam, mu);
        std::cout << "dP/dF analytical:\n";
        std::cout << ToString(dP_dF_mult_trick) << std::endl << std::endl;

        TensorDiffStats(dP_dF_FD, dP_dF_mult_trick);
    }
}


void Optimization()
{
    static int layer = 0;

    if (ImGui::TreeNode("Optimization Setup"))
    {
        if (ImGui::Button("Load Optimization Input from JSON"))
        {
            std::ifstream ifs;
            ifs.open("optimization_input_test.json");

            if (ifs.good()) {
                cereal::JSONInputArchive iarchive(ifs); // Create an input archive

                iarchive(opt_input); // Read the data from the archive
                ifs.close();
            }
            else {
                std::cout << "couldn't open input json" << std::endl;
            }
        }

        if (ImGui::Button("Save Optimization Input from JSON"))
        {
            std::ofstream ofs;
            ofs.open("optimization_input_test.json");

            if (ofs.good()) {
                cereal::JSONOutputArchive oarchive(ofs); // Create an output archive

                oarchive(opt_input); // Write the data to the archive
            } // archive goes out of scope, ensuring all contents are flushed
            ofs.close();
        }

        opt_input.ImGui();

        if (ImGui::Button("Construct computation graph from optimization input"))
        {
            LoadCompGraph(opt_input, comp_graph, &ps_point_cloud, &ps_target_point_cloud, &ps_grid);
        }
        ImGui::TreePop();
    }


    if (!comp_graph)
        ImGui::BeginDisabled();


    // Might need to do some checks to make sure stuff is loaded before button is pressed
    if (ImGui::TreeNode("Advanced Visualization"))
    {
        if (ImGui::Button("Add Target Grid Masses to Grid Visualization"))
        {
            auto target_grid_masses = comp_graph->target_grid->GetNodeMasses();
            ps_grid->addScalarQuantity("target grid masses", target_grid_masses);

        }
        if (ImGui::Button("Add control point elastic energies to visualization"))
        {
            auto elastic_energies = comp_graph->layers[layer].point_cloud->GetPointElasticEnergies();
            ps_point_cloud->addScalarQuantity("elastic energies", elastic_energies);
        }

        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Test Buttons"))
    {
        if (ImGui::Button("Test binary point cloud file writing"))
        {
            comp_graph->layers[layer].point_cloud->WriteEntirePointCloudToBinaryFile("test_binary.mpm");
            comp_graph->layers[layer].point_cloud->WriteEntirePointCloudToFile("test_text.mpm");
        }
        ImGui::TreePop();
    }


    static int point_index = 0;
    ImGui::InputInt("point index", &point_index);
    if (ImGui::Button("Remove point")) {
        comp_graph->layers.begin()->point_cloud->RemovePoint(point_index);
        auto point_positions = comp_graph->layers.begin()->point_cloud->GetPointPositions();

        ps_point_cloud->updatePointPositions(point_positions);
        ps_point_cloud->refresh();
    }


    if (ImGui::Button("Set initial deformation gradients"))
    {

        auto mpm_pc = comp_graph->layers.front().point_cloud;
        for (size_t i = 0; i < mpm_pc->points.size(); i++) {
            // ROTATION MATRIX
            double a = 45.0 * 0.0174533;
            Mat3 F;
            F << 1, 0, 0, 0, cos(a), -sin(a), 0, sin(a), cos(a);

            //Mat3 F;
            //F.setIdentity();
            F *= 0.8;
            mpm_pc->points[i].F = F;
        }
    }

    static double young_mod = 400.0;
    ImGui::InputDouble("Young's Modulus", &young_mod);
    static double poisson = 0.49;
    ImGui::InputDouble("Poisson's Ratio", &poisson);
    if (ImGui::Button("Calculate Lame Parameters"))
    {
        double lam, mu;
        CalculateLameParameters(young_mod, poisson, lam, mu);
        std::cout << "lam = " << lam << ", mu = " << mu << std::endl;
        auto mpm_pc = comp_graph->layers.front().point_cloud;
        for (size_t i = 0; i < mpm_pc->points.size(); i++) {
            mpm_pc->points[i].lam = lam;
            mpm_pc->points[i].mu = mu;
        }
    }

    if (ImGui::Button("Finite Differences test dLdF"))
    {
        comp_graph->FiniteDifferencesGradientTest(opt_input.num_timesteps, 0);
    }

    if (ImGui::Button("Optimize Control Sequence"))
    {
        auto begin_clock = std::chrono::steady_clock::now();
        for (int i = 0; i < opt_input.num_animations; i++)
        {
            auto curr_begin_clock = std::chrono::steady_clock::now();
            std::cout << "**********OPTIMIZING ANIMATION INTERVAL: " << i << "************" << std::endl;
            comp_graph->layers.front().point_cloud = comp_graph->layers.back().point_cloud;
            comp_graph->layers.resize(1);


            comp_graph->OptimizeDefGradControlSequence(
                opt_input.num_timesteps,
                opt_input.dt,
                opt_input.drag,
                opt_input.f_ext,
                opt_input.control_stride,
                opt_input.max_gd_iters,
                opt_input.max_ls_iters,
                opt_input.initial_alpha,
                opt_input.gd_tol
            );

            // RENDER
            for (size_t t = 0; t < comp_graph->layers.size(); t++)
            {
                std::string mpm_output_folder = "output\\";
                std::string number_str = LeadingZerosNumberStr(i * opt_input.num_timesteps + t, 6);

                auto point_positions = comp_graph->layers[t].point_cloud->GetPointPositions();
                ps_point_cloud->updatePointPositions(point_positions);

                // VISUALIZE ELASTIC ENERGIES
                auto elastic_energies = comp_graph->layers[t].point_cloud->GetPointElasticEnergies();
                ps_point_cloud->addScalarQuantity("elastic energies", elastic_energies);


                std::string png_output_path = mpm_output_folder + "screenshot_" + number_str + ".png";
                polyscope::screenshot(png_output_path, false);

                // SAVE POINT DATA TO FILE
                std::string mpm_output_path = mpm_output_folder + "mpm_points_" + number_str + ".mpm";
                comp_graph->layers[t].point_cloud->WriteEntirePointCloudToFile(mpm_output_path);
            }

            auto curr_end_clock = std::chrono::steady_clock::now();
            std::cout << "Animation interval " << i << " took " << std::chrono::duration_cast<std::chrono::seconds>(curr_end_clock - curr_begin_clock).count() << " seconds." << std::endl;
            std::cout << "Full animation took " << std::chrono::duration_cast<std::chrono::seconds>(curr_end_clock - begin_clock).count() << " seconds so far." << std::endl;

        }
    }




    if (ImGui::InputInt("layer", &layer) && comp_graph) {
        if (layer >= (int)comp_graph->layers.size())
            layer = (int)comp_graph->layers.size() - 1;
        if (layer < 0)
            layer = 0;

        auto point_positions = comp_graph->layers[layer].point_cloud->GetPointPositions();
        ps_point_cloud->updatePointPositions(point_positions);
    }




    if (!comp_graph)
        ImGui::EndDisabled();
}

void menuCallback()
{
    ImGui::PushItemWidth(100);

    if (ImGui::TreeNode("Gradient Correctness Checking"))
    {
        CheckGradients();
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("RT MPM GUI"))
    {
        RealTimeMPMImGUI();
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Optimization"))
    {
        Optimization();
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Visualization/Animation"))
    {
        ImGui::PushItemWidth(300);

        if (ImGui::Button("Visualize Grid"))
        {
            // MPM Grid
            std::cout << "generating mpm grid..." << std::endl;
            int grid_dims[3];
            for (int i = 0; i < 3; i++) {
                grid_dims[i] = (int)std::ceil((opt_input.grid_max_point[0] - opt_input.grid_min_point[0]) / opt_input.grid_dx);
            }
            mpm_grid = std::make_shared<Grid>(grid_dims[0], grid_dims[1], grid_dims[2], opt_input.grid_dx, opt_input.grid_min_point);
            // grid nodes
            auto grid_points = mpm_grid->GetNodePositions();
            ps_grid = polyscope::registerPointCloud(PS_SIM_GRID, grid_points);
            ps_grid->setPointRadius(mpm_grid->dx / 500.0);
            ps_grid->setPointRenderMode(polyscope::PointRenderMode::Sphere);
        }


        static std::string pc_data_folder = "experiments/big_sca_demo/output_bob_to_spot/"; //mpm_points_000000.mpm
        ImGui::InputText("Point cloud data folder", &pc_data_folder);

        static std::string pc_file = "mpm_points_000000.mpm";
        ImGui::InputText("Point cloud file", &pc_file);

        if (ImGui::Button("Load point cloud"))
        {
            mpm_point_cloud = std::make_shared<PointCloud>();
            mpm_point_cloud->ReadEntirePointCloudFromFile(pc_data_folder + pc_file);
            auto positions = mpm_point_cloud->GetPointPositions();
            ps_point_cloud = polyscope::registerPointCloud(PS_POINT_CLOUD_1, positions);
            double point_dx = opt_input.grid_dx / (double)opt_input.points_per_cell_cuberoot;
            ps_point_cloud->setPointRadius(point_dx / 50.0);
            ps_point_cloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
        }

        ImGui::PopItemWidth();
        ImGui::TreePop();
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
    polyscope::state::userCallback = menuCallback;
    //polyscope::state::userCallback = realtimeCallback;

    // Pass control flow to polyscope, displaying the interactive window.
    // Function will return when user closes the window.
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::show();
}