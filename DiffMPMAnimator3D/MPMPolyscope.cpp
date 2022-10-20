#include "MPMPolyscope.h"
#include "igl/readOBJ.h"
#include "GeometryLoading.h"
#include "ForwardSimulation.h"

bool LoadMPMPointCloudFromObj(
    std::string obj_path,
    std::shared_ptr<PointCloud>& mpm_point_cloud,
    double point_dx,
    double density,
    double lam,
    double mu
)
{
    std::cout << "reading " << obj_path << "..." << std::endl;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (!igl::readOBJ(obj_path, V, F))
        return false;
    Vec3 min_point = Vec3(DBL_MAX, DBL_MAX, DBL_MAX);
    Vec3 max_point = Vec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    for (size_t i = 0; i < (size_t)V.rows(); i++) {
        for (size_t j = 0; j < 3; j++) {
            min_point[j] = std::fmin(min_point[j], V(i, j));
            max_point[j] = std::fmax(max_point[j], V(i, j));
        }
    }
    std::cout << obj_path << "'s max point is: " << max_point.transpose() << std::endl;
    std::cout << obj_path << "'s min point is: " << min_point.transpose() << std::endl;

    std::cout << "generating point cloud..." << std::endl;
    std::vector<Vec3> points = GeometryLoading::GeneratePointCloudFromWatertightTriangleMesh(V, F, min_point, max_point, point_dx);


    // MPM Point Cloud
    std::cout << "generating mpm points..." << std::endl;
    mpm_point_cloud = std::make_shared<PointCloud>();
    mpm_point_cloud->points.resize(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        mpm_point_cloud->points[i].x = points[i];
        mpm_point_cloud->points[i].v = Vec3(0.0, 0.0, 0.0);
        mpm_point_cloud->points[i].F = Mat3::Identity();
        mpm_point_cloud->points[i].m = point_dx * point_dx * point_dx * density;
        mpm_point_cloud->points[i].lam = lam;
        mpm_point_cloud->points[i].mu = mu;
    }

    return true;
}

bool LoadScene(const SceneInput& scene_input,
    std::shared_ptr<PointCloud>& mpm_point_cloud, std::shared_ptr<Grid>& mpm_grid, 
    polyscope::PointCloud** polyscope_point_cloud,
    polyscope::PointCloud** polyscope_grid)
{
    double point_dx = scene_input.grid_dx / (double)scene_input.points_per_cell_cuberoot;

    
    if (!LoadMPMPointCloudFromObj(scene_input.mpm_input_mesh_path, mpm_point_cloud, point_dx, scene_input.p_density, scene_input.lam, scene_input.mu))
        return false;
    


    // MPM Grid
    std::cout << "generating mpm grid..." << std::endl;
    int grid_dims[3];
    for (int i = 0; i < 3; i++) {
        grid_dims[i] = (int)std::ceil((scene_input.grid_max_point[0] - scene_input.grid_min_point[0]) / scene_input.grid_dx);
    } 
    mpm_grid = std::make_shared<Grid>(grid_dims[0], grid_dims[1], grid_dims[2], scene_input.grid_dx, scene_input.grid_min_point);


    // calculate volume
    SingleThreadMPM::CalculatePointCloudVolumes(*mpm_point_cloud, *mpm_grid);

    // Polyscope point cloud
    std::cout << "registering point cloud..." << std::endl;
    auto points = mpm_point_cloud->GetPointPositions();
    *polyscope_point_cloud = polyscope::registerPointCloud(scene_input.mpm_input_mesh_path, points);
    (*polyscope_point_cloud)->setPointRadius(point_dx / 50.0); // Not sure how this radius is actually calculated, this is hand tuned
    (*polyscope_point_cloud)->setPointRenderMode(polyscope::PointRenderMode::Sphere);

    auto grid_points = mpm_grid->GetNodePositions();
    *polyscope_grid = polyscope::registerPointCloud("mpm_grid", grid_points);
    (*polyscope_grid)->setPointRadius(mpm_grid->dx / 500.0);
    (*polyscope_grid)->setPointRenderMode(polyscope::PointRenderMode::Sphere);

	return true;
}

bool LoadCompGraph(
    const SceneInput& scene_input,
    std::shared_ptr<CompGraph>& comp_graph,
    polyscope::PointCloud** polyscope_input_point_cloud,
    polyscope::PointCloud** polyscope_target_point_cloud,
    polyscope::PointCloud** polyscope_grid
)
{
    std::shared_ptr<PointCloud> initial_point_cloud;
    std::shared_ptr<Grid> grid;


    // Load simulation scene
    if (!LoadScene(scene_input, initial_point_cloud, grid, polyscope_input_point_cloud, polyscope_grid))
        return false;

    
    // also load the target point cloud
    double point_dx = scene_input.grid_dx / (double)scene_input.points_per_cell_cuberoot;
    std::shared_ptr<PointCloud> target_point_cloud;
    if (!LoadMPMPointCloudFromObj(scene_input.mpm_target_mesh_path, target_point_cloud, point_dx, scene_input.p_density, scene_input.lam, scene_input.mu))
        return false;

    // calculate volume
    SingleThreadMPM::CalculatePointCloudVolumes(*target_point_cloud, *grid);

    // Polyscope point cloud
    std::cout << "registering point cloud..." << std::endl;
    auto points = target_point_cloud->GetPointPositions();
    *polyscope_target_point_cloud = polyscope::registerPointCloud(scene_input.mpm_target_mesh_path, points);
    (*polyscope_target_point_cloud)->setPointRadius(point_dx / 50.0); // Not sure how this radius is actually calculated, this is hand tuned
    (*polyscope_target_point_cloud)->setPointRenderMode(polyscope::PointRenderMode::Sphere);

    // Then turn it into a target grid
    auto target_grid = std::make_shared<Grid>(*grid);
    SingleThreadMPM::P2G(*target_point_cloud, *target_grid, 0.0, 0.0);


    // Finally make the comp graph
    comp_graph = std::make_shared<CompGraph>(initial_point_cloud, grid, target_grid);

    return true;
}