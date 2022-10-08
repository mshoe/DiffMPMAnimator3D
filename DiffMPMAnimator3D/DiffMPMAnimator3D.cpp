// DiffMPMAnimator3D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"



// C++ includes
#include <iostream>
#include <sstream>
using namespace std;

// DiffMPMLib3D
#include "MPMPolyscope.h"
#include "ForwardSimulation.h"

SceneInput scene_input;

// FOR REAL TIME
std::shared_ptr<PointCloud> mpm_point_cloud = nullptr;
std::shared_ptr<Grid> mpm_grid = nullptr;
polyscope::PointCloud* ps_point_cloud = nullptr;
polyscope::PointCloud* ps_grid = nullptr;

void realtimeCallback() 
{
    ImGui::PushItemWidth(100);

    double screenshot_interval = 1.0 / 30.0;

    static double time = 0.0;
    static double time_until_next_screenshot = 0.0;
    static int frame = 0;

    static bool scene_loaded = false;
    if (scene_loaded)
        ImGui::BeginDisabled();
    if (ImGui::Button("Load Scene"))
    {
        LoadScene(scene_input, mpm_point_cloud, mpm_grid, &ps_point_cloud, &ps_grid);
        scene_loaded = true;
    }
    if (scene_loaded)
        ImGui::EndDisabled();

    if (!scene_loaded)
    {
        ImGui::PopItemWidth();
        return;
    }

    if (ImGui::Button("timestep")) {
        // executes when button is pressed
        ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, scene_input.dt, scene_input.drag, scene_input.f_ext);
        time += scene_input.dt;
        frame++;

        auto point_positions = mpm_point_cloud->GetPointPositions();
        ps_point_cloud->updatePointPositions(point_positions);
    }

    if (ImGui::Button("screenshot")) {
        polyscope::screenshot();
    }

    
    static bool playing = false;
    ImGui::Checkbox("playing", &playing);
    if (playing) {
        ForwardTimeStep(*mpm_point_cloud, *mpm_point_cloud, *mpm_grid, scene_input.dt, scene_input.drag, scene_input.f_ext);
        time += scene_input.dt;
        time_until_next_screenshot -= scene_input.dt;
        frame++;

        


        
        if (time_until_next_screenshot < 0.0) {
            auto point_positions = mpm_point_cloud->GetPointPositions();
            ps_point_cloud->updatePointPositions(point_positions);
            time_until_next_screenshot = screenshot_interval;
            polyscope::screenshot();
        }

        if (time > 10.0) {
            playing = false;
        }
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


    ImGui::PopItemWidth();
}



// FOR OPTIMIZATION
std::shared_ptr<CompGraph> comp_graph = nullptr;

void optimizationCallback()
{
    ImGui::PushItemWidth(100);

    if (ImGui::Button("Optimize Control Sequence"))
    {
        LoadCompGraph(scene_input, comp_graph, &ps_point_cloud, &ps_grid);

        comp_graph->OptimizeDefGradControlSequence(
            60,
            scene_input.dt,
            scene_input.drag,
            scene_input.f_ext,
            10,
            25,
            15,
            0.025,
            1e-3
        );
    }


    static int layer = 0;
    if (ImGui::InputInt("layer", &layer) && comp_graph) {
        if (layer >= (int)comp_graph->layers.size())
            layer = (int)comp_graph->layers.size() - 1;
        if (layer < 0)
            layer = 0;
        
        auto point_positions = comp_graph->layers[layer].point_cloud->GetPointPositions();
        ps_point_cloud->updatePointPositions(point_positions);
    }

    ImGui::PopItemWidth();
}


int main()
{
    std::cout << "Hello World!\n";
    

    // Initialize polyscope, creating graphics contexts and constructing a window.
    // Should be called exactly once.
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::init();

    /*
    * build visualizations, here or in distant code
    *
    */
    

    // Specify the callback
    polyscope::state::userCallback = optimizationCallback;

    // Pass control flow to polyscope, displaying the interactive window.
    // Function will return when user closes the window.
    polyscope::show();

    // More of your code
    // ...

    

    // Show again. Data is preserved between calls to show()
    // unless explicitly removed.
    //polyscope::show();

    //todo: make mpm library
}