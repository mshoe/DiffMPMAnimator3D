#include "pch.h"
#include "Elasticity.h"
#include "qrsvd/ImplicitQRSVD.h"
#include <unsupported/Eigen/KroneckerProduct>

Tensor3x3x3x3 d_JFit_dF_FD(const Mat3& F)
{
    // return d (J*Fit)_ab / d F_ij

    Tensor3x3x3x3 ret;

    auto CalcJFit = [](const Mat3& _F)
    {
        double J = _F.determinant();
        Mat3 Fit = _F.inverse().transpose();
        Mat3 ret = J * Fit;
        return ret;
    };
    /*std::cout << "F\n" << F << std::endl;
    std::cout << "inv(F)\n" << F.inverse() << std::endl;
    std::cout << "det(F)\n" << F.determinant() << std::endl;
    std::cout << "inv(F)^T\n" << F.inverse().transpose() << std::endl;
    std::cout << "JFit inline\n" << F.determinant() * F.inverse().transpose() << std::endl;
    double J = F.determinant();
    Mat3 Fit = F.inverse().transpose();
    auto JFit = J * Fit;
    std::cout << "JFit inline auto\n" << JFit << std::endl;
    std::cout << "JFit lambda\n" << CalcJFit(F) << std::endl;*/
    
    Mat3 temp_F = F;
    double delta = 1e-9;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double originalValue = temp_F(i, j);

            temp_F(i, j) = originalValue + delta;
            Mat3 JFit_forward = CalcJFit(temp_F);
            //std::cout << "JFit_forward\n" << JFit_forward << std::endl;

            temp_F(i, j) = originalValue - delta;
            Mat3 JFit_backward = CalcJFit(temp_F);
            //std::cout << "JFit_backward\n" << JFit_backward << std::endl;

            temp_F(i, j) = originalValue;
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++)
                {
                    double l1 = JFit_forward(a, b);
                    double l2 = JFit_backward(a, b);
                    ret[a][b](i, j) = (l1 - l2) / (2.0 * delta);
                }
            }
        }
    }

    return ret;
}

Mat3 PK_FixedCorotatedElasticity(const Mat3& F, double lam, double mu)
{
    Mat3 R, S;
    JIXIE::polarDecomposition(F, R, S);

    double J = F.determinant();

    return 2.0 * mu * (F - R) + lam * (J - 1.0) * J * F.inverse().transpose();
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
    Mat3 Finv = F.inverse();
    Mat3 Fit = Finv.transpose();

    Mat3 ret = Mat3::Zero();

    // 2 mu dF 
    ret += 2.0 * mu * dF;

    // + lam J F^-T (J F^-T : dF)
    ret += lam * J * Fit * InnerProduct(J * Fit, dF);
    
    // + lam (J - 1) d(J F^-T)
    // d(J F^-T) derivation in my notes
    //Mat3 cofactor_F = F.adjoint().transpose();
    //Mat3 E = Mat3::Zero(); // unit matrix
    //Mat3 term = Mat3::Zero();
    //for (int i = 0; i < 3; i++) {
    //    for (int j = 0; j < 3; j++) {
    //        E.setZero();
    //        E(i, j) = 1.0;
    //        Mat3 temp_right = J * (-Finv * E * Finv).transpose();
    //        Mat3 temp_left = cofactor_F(i, j) * Fit;
    //        Mat3 temp = temp_left + temp_right;
    //        term(i, j) = InnerProduct(temp, dF);
    //    }
    //}
    //ret += lam * (J - 1.0) * term;


    // + lam (J - 1) d(J F^-T)
    // d(J F^-T) Derivation done in sympy
    //Eigen::KroneckerProduct<Mat3, Mat3> dJFit(J * Fit, -Finv);
    ////dJFit.
    //std::cout << dJFit << std::endl << std::endl;
    /*double F00, F01, F02, F10, F11, F12, F20, F21, F22;
    F00 = F(0, 0);
    F01 = F(0, 1);
    F02 = F(0, 2);
    F10 = F(1, 0);
    F11 = F(1, 1);
    F12 = F(1, 2);
    F20 = F(2, 0);
    F21 = F(2, 1);
    F22 = F(2, 2);
    std::vector<std::vector<Mat3>> dJFit_dF(3, std::vector<Mat3>(3, Mat3::Identity()));
    
    
    Mat3 temp_00;
    temp_00 << 0, 0, 0, 0, F22, -F12, 0, -F12, F11;
    dJFit_dF[0][0] = temp_00;
    
    Mat3 temp_01;
    temp_01 << 0, -F22, F12, -F22, 0, F02, F12, F02, -2 * F01;
    dJFit_dF[0][1] = temp_01;
    
    Mat3 temp_02;
    temp_02 << 0, F12, -F11, F12, -2 * F02, F01, -F11, F01, 0;
    dJFit_dF[0][2] = temp_02;

    Mat3 temp_10;
    temp_10 << 0, -F22, F12, -F22, 0, F02, F12, F02, -2 * F01;
    dJFit_dF[1][0] = temp_10;
    
    Mat3 temp_11;
    temp_11 << F22, 0, -F02, 0, 0, 0, -F02, 0, F00;
    dJFit_dF[1][1] = temp_11;

    Mat3 temp_12;
    temp_12 << -2 * F12, F02, F01, F02, 0, -F00, F01, -F00, 0;
    dJFit_dF[1][2] = temp_12;

    Mat3 temp_20;
    temp_20 << 0, F12, -F11, F12, -2 * F02, F01, -F11, F01, 0;
    dJFit_dF[2][0] = temp_20;

    Mat3 temp_21;
    temp_21 << -2 * F12, F02, F01, F02, 0, -F00, F01, -F00, 0;
    dJFit_dF[2][1] = temp_21;

    Mat3 temp_22;
    temp_22 << F11, -F01, 0, -F01, F00, 0, 0, 0, 0;
    dJFit_dF[2][2] = temp_22;*/

    Tensor3x3x3x3 dJFit_dF = d_JFit_dF_FD(F);

    Mat3 term = Mat3::Zero();
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) 
        {
            term(a, b) = 0.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    term(a, b) += dJFit_dF[a][b](i, j) * dF(i, j);
                }
            }
        }
    }
    ret += lam * (J - 1.0) * term;

    // Analytical
    //ret += lam * (J - 1.0) * dF.adjoint().transpose();

    // Now time to calculate dR
    // Derivation done in sympy
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


Tensor3x3x3x3 d2_FCE_psi_dF2_FD(const Mat3& F, double lam, double mu)
{
    // Return tensor of indices: a, b, i, j,  for dP[a][b] / dF[i][j]

    Tensor3x3x3x3 ret;

    double delta = 1e-9;
    Mat3 temp_F = F;
    Mat3 P = PK_FixedCorotatedElasticity(temp_F, lam, mu);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double originalValue = temp_F(i, j);

            temp_F(i, j) = originalValue + delta;
            Mat3 P_forward = PK_FixedCorotatedElasticity(temp_F, lam, mu);

            temp_F(i, j) = originalValue - delta;
            Mat3 P_backward = PK_FixedCorotatedElasticity(temp_F, lam, mu);

            temp_F(i, j) = originalValue;
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) 
                {
                    double l1 = P_forward(a, b);
                    double l2 = P_backward(a, b);
                    ret[a][b](i, j) = (l1 - l2) / (2.0 * delta);
                }
            }
        }
    }

    return ret;
}

Tensor3x3x3x3 d2_FCE_psi_dF2_mult_trick(const Mat3& F, double lam, double mu)
{
    Tensor3x3x3x3 ret;
    Mat3 E;// = Mat3::Identity();
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            E.setZero();
            E(a, b) = 1.0;
            ret[a][b] = d2_FCE_psi_dF2_mult_by_dF(F, lam, mu, E);
        }
    }
    return ret;
}
