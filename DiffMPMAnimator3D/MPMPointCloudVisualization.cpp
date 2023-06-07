#include "MPMPointCloudVisualization.h"
//#include "MPMPolyscope.h"
#include "ForwardSimulation.h"
#include "Interpolation.h"
#include "polyscope/scalar_quantity.h"

std::vector<double> GetPointCloudLocalMassFields(const DiffMPMLib3D::PointCloud& mpm_point_cloud, DiffMPMLib3D::Grid& grid)
{
	using namespace DiffMPMLib3D;
    // project mass to grid
	grid.ResetValues();
    SingleThreadMPM::P2G(mpm_point_cloud, grid, 0.0, 0.0);

    std::vector<double> local_mass_field;
	local_mass_field.resize(mpm_point_cloud.points.size());

    for (size_t p = 0; p < mpm_point_cloud.points.size(); p++) {
		const MaterialPoint& mp = mpm_point_cloud.points[p];
		local_mass_field[p] = 0.0;

        // G2P to get local mass field of particles
		double dx = grid.dx;

		Vec3 relative_point = mp.x - grid.min_point;
		int bot_left_index[3];
		for (size_t i = 0; i < 3; i++) {
			bot_left_index[i] = (int)std::floor(relative_point[i] / grid.dx) - 1;
		}
		for (int i = 0; i <= 3; i++) {
			for (int j = 0; j <= 3; j++) {
				for (int k = 0; k <= 3; k++) {
					int index[3] = {
						bot_left_index[0] + i,
						bot_left_index[1] + j,
						bot_left_index[2] + k
					};
					if (0 <= index[0] && index[0] < grid.dim_x &&
						0 <= index[1] && index[1] < grid.dim_y &&
						0 <= index[2] && index[2] < grid.dim_z)
					{
						const GridNode& node = grid.nodes[index[0]][index[1]][index[2]];
						Vec3 xg = node.x;
						Vec3 xp = mp.x;
						Vec3 dgp = xg - xp;
						double wgp = CubicBSpline(dgp[0] / grid.dx) * CubicBSpline(dgp[1] / grid.dx) * CubicBSpline(dgp[2] / grid.dx);

						local_mass_field[p] += wgp * node.m;
					}
				}
			}
		}

		// post process a bit
		if (local_mass_field[p] < 1.0) {
			local_mass_field[p] = local_mass_field[p] + 100.0 * std::log(local_mass_field[p]);
		}
    }

    return local_mass_field;
}

void UpdatePolyscopePointCloudProperties(polyscope::PointCloud** ps_point_cloud, std::shared_ptr<const DiffMPMLib3D::PointCloud> mpm_point_cloud)
{
    auto positions = mpm_point_cloud->GetPointPositions();
    (*ps_point_cloud)->updatePointPositions(positions);

    auto elastic_energies = mpm_point_cloud->GetPointElasticEnergies();
    (*ps_point_cloud)->addScalarQuantity("elastic_energies", elastic_energies);

    auto masses = mpm_point_cloud->GetPointMasses();
    (*ps_point_cloud)->addScalarQuantity("m", masses);

    auto volumes = mpm_point_cloud->GetPointVolumes();
    (*ps_point_cloud)->addScalarQuantity("vol", volumes);

    auto determinants = mpm_point_cloud->GetPointDeterminants();
    (*ps_point_cloud)->addScalarQuantity("J", determinants);
}

void UpdatePolyscopePointCloudMassField(polyscope::PointCloud** ps_point_cloud, std::shared_ptr<const DiffMPMLib3D::PointCloud> mpm_point_cloud, std::shared_ptr<DiffMPMLib3D::Grid> grid, double min_val, double max_val)
{
	auto positions = mpm_point_cloud->GetPointPositions();
	(*ps_point_cloud)->updatePointPositions(positions);

	auto local_mass_field = GetPointCloudLocalMassFields(*mpm_point_cloud, *grid);
	auto q = (*ps_point_cloud)->addScalarQuantity("local mass field", local_mass_field);
	q->setMapRange({ min_val, max_val });
	//(*ps_point_cloud)->setMap
}