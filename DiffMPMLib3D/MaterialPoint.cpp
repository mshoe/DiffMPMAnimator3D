#include "pch.h"
#include "MaterialPoint.h"

void DiffMPMLib3D::MaterialPoint::ResetGradients()
{
	dLdx.setZero();
	dLdv.setZero();
	dLdF.setZero();
	dLdC.setZero();
	dLdP.setZero();

	dLdv_next.setZero();
	dLdC_next.setZero();

	dLdm = 0.0;
	dLdvol = 0.0;
}

void DiffMPMLib3D::MaterialPoint::WriteEntirePointToFile(std::ofstream& ofs)
{
    ofs << x(0) << " " << x(1) << " " << x(2) << std::endl;
    ofs << v(0) << " " << v(1) << " " << v(2) << std::endl;
    ofs << F(0, 0) << " " << F(0, 1) << " " << F(0, 2) << " " << F(1, 0) << " " << F(1, 1) << " " << F(1, 2) << " " << F(2, 0) << " " << F(2, 1) << " " << F(2, 2) << std::endl;
    ofs << C(0, 0) << " " << C(0, 1) << " " << C(0, 2) << " " << C(1, 0) << " " << C(1, 1) << " " << C(1, 2) << " " << C(2, 0) << " " << C(2, 1) << " " << C(2, 2) << std::endl;
    ofs << P(0, 0) << " " << P(0, 1) << " " << P(0, 2) << " " << P(1, 0) << " " << P(1, 1) << " " << P(1, 2) << " " << P(2, 0) << " " << P(2, 1) << " " << P(2, 2) << std::endl;
    ofs << m << std::endl;
    ofs << vol << std::endl;
    ofs << lam << std::endl;
    ofs << mu << std::endl;
    ofs << dFc(0, 0) << " " << dFc(0, 1) << " " << dFc(0, 2) << " " << dFc(1, 0) << " " << dFc(1, 1) << " " << dFc(1, 2) << " " << dFc(2, 0) << " " << dFc(2, 1) << " " << dFc(2, 2) << std::endl;
    
    // The gradients do not need to be stored in order to load the point cloud,
    // but we will store dLdx, dLdv, and dLdF for visualization purposes
    ofs << dLdx(0) << " " << dLdx(1) << " " << dLdx(2) << std::endl;
    ofs << dLdv(0) << " " << dLdv(1) << " " << dLdv(2) << std::endl;
    ofs << dLdF(0, 0) << " " << dLdF(0, 1) << " " << dLdF(0, 2) << " " << dLdF(1, 0) << " " << dLdF(1, 1) << " " << dLdF(1, 2) << " " << dLdF(2, 0) << " " << dLdF(2, 1) << " " << dLdF(2, 2) << std::endl;
}

void DiffMPMLib3D::MaterialPoint::ReadEntirePointFromFile(std::ifstream& ifs)
{
    ifs >> x(0) >> x(1) >> x(2);
    ifs >> v(0) >> v(1) >> v(2);
    ifs >> F(0, 0) >> F(0, 1) >> F(0, 2) >> F(1, 0) >> F(1, 1) >> F(1, 2) >> F(2, 0) >> F(2, 1) >> F(2, 2);
    ifs >> C(0, 0) >> C(0, 1) >> C(0, 2) >> C(1, 0) >> C(1, 1) >> C(1, 2) >> C(2, 0) >> C(2, 1) >> C(2, 2);
    ifs >> P(0, 0) >> P(0, 1) >> P(0, 2) >> P(1, 0) >> P(1, 1) >> P(1, 2) >> P(2, 0) >> P(2, 1) >> P(2, 2);
    ifs >> m;
    ifs >> vol;
    ifs >> lam;
    ifs >> mu;
    ifs >> dFc(0, 0) >> dFc(0, 1) >> dFc(0, 2) >> dFc(1, 0) >> dFc(1, 1) >> dFc(1, 2) >> dFc(2, 0) >> dFc(2, 1) >> dFc(2, 2);
    ifs >> dLdx(0) >> dLdx(1) >> dLdx(2);
    ifs >> dLdv(0) >> dLdv(1) >> dLdv(2);
    ifs >> dLdF(0, 0) >> dLdF(0, 1) >> dLdF(0, 2) >> dLdF(1, 0) >> dLdF(1, 1) >> dLdF(1, 2) >> dLdF(2, 0) >> dLdF(2, 1) >> dLdF(2, 2);
}
