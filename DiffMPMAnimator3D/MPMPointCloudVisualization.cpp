#include "MPMPointCloudVisualization.h"
//#include "MPMPolyscope.h"

void UpdatePolyscopePointCloudProperties(polyscope::PointCloud** ps_point_cloud, std::shared_ptr<DiffMPMLib3D::PointCloud> mpm_point_cloud)
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
