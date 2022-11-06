#include "pch.h"
#include "GeometryLoading.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/signed_distance.h"

std::vector<DiffMPMLib3D::Vec3> DiffMPMLib3D::GeometryLoading::GeneratePointCloudFromWatertightTriangleMesh(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi& F,
    Vec3 min_point,
    Vec3 max_point,
    double sampling_dx)
{
    using namespace Eigen;

    // 1. generate our uniform grid sample points
    
    // determine dimensions of grid sampling
    int dims[3];
    for (int i = 0; i < 3; i++) {
        dims[i] = std::ceil((max_point[i] - min_point[i]) / sampling_dx);
    }

    // construct sample points
    MatrixXd P;
    P.resize(dims[0] * dims[1] * dims[2], 3);
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            for (int k = 0; k < dims[2]; k++) {
                P.row(i * dims[1] * dims[2] + j * dims[2] + k) = min_point + sampling_dx * Vec3(double(i), double(j), double(k));
            }
        }
    }


    // 2. get squared distance queries from sample points to mesh
    VectorXd S;
    VectorXi I;
    MatrixXd C;
    MatrixXd N;

    // 3. get signed distances
    // Choose type of signing to use
    igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
    igl::signed_distance(P, V, F, sign_type, S, I, C, N);


    // 4. store all points with negative signed distance
    std::vector<Vec3> points;
    for (int i = 0; i < P.rows(); i++) {
        if (S[i] <= 0.0) {
            points.push_back(P.row(i));
        }
    }

    return points;
}
