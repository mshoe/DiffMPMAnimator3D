#include "pch.h"
#include "Elasticity.h"
#include "qrsvd/ImplicitQRSVD.h"

Mat3 PK_FixedCorotatedElasticity(const Mat3& F, double lam, double mu)
{
    Mat3 R, S;
    JIXIE::polarDecomposition(F, R, S);

    double J = F.determinant();

    return 2.0 * mu * (F - R) + lam * (J - 1.0) * J * F.inverse().transpose();
}
