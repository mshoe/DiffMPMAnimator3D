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

// calculates the derivative of X.inverse().transpose() wrt. X(i, j)
// Given Xinv = X.inverse()
Mat3 d_Xit_dX_ij(const Mat3& Xinv, int i, int j)
{
    assert(0 <= i && i < 3);
    assert(0 <= j && j < 3);

    Mat3 ret;

    for (int l = 0; l < 3; l++) {
        for (int k = 0; k < 3; k++) {
            ret(l, k) = Xinv(k, i)* Xinv(j, l);
        }
    }

    return ret;
}


Mat3 d2_FCE_psi_dF2_mult_by_dF(const Mat3& F, double lam, double mu, const Mat3& dF)
{
    /*
    * Polar decomp: F = R S
    * 
    * (d2psi / dFdF)(F) : dF  =
    * 
    *   2 mu dF 
    * - 2 mu dR
    * + lam J F^-T (J F^-T : dF)
    * + lam (J - 1) d(J F^-T)
    * 
    */


    Mat3 R, S;
    JIXIE::polarDecomposition(F, R, S);
    double J = F.determinant();
    Mat3 Fit = F.inverse().transpose();

    Mat3 ret = Mat3::Zero();

    // 2 mu dF 
    ret += 2.0 * mu * dF;

    // + lam J F^-T (J F^-T : dF)
    ret += lam * J * Fit * InnerProduct(J * Fit, dF);
    
    //ret += lam * (J - 1.0) * CofactorMatrix(dF);
    


    // Now time to calculate dR
    double s00, s01, s02, s11, s12, s22;
    s00 = S(0, 0);
    s01 = S(0, 1);
    s02 = S(0, 2);
    s11 = S(1, 1);
    s12 = S(1, 2);
    s22 = S(2, 2);
    Mat3 A;
    A.row(0) = Vec3(s02         , s12       , -(s00 + s11)  );
    A.row(1) = Vec3(-s01        , s00 + s22 , -s12          );
    A.row(2) = Vec3(-(s11 + s22), s01       , s02           );

    Mat3 temp = R.transpose() * dF - dF.transpose() * R;
    Vec3 b = Vec3(temp(0, 1), temp(0, 2), temp(1, 2));

    // Let x = R.T * dR
    // solve Ax = b for x
    Vec3 x = A.colPivHouseholderQr().solve(b);

    // dR = R(R.T * dR)
    Mat3 dR;
    dR.row(0) = Vec3(0      , -x(2) , x(1)  );
    dR.row(1) = Vec3(x(2)   , 0     , -x(0) );
    dR.row(2) = Vec3(-x(1)  , x(0)  , 0.0   );
    dR = R * dR;

    ret -= 2.0 * mu * dR;

    return ret;
}
