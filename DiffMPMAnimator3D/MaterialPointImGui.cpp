#include "MaterialPointImGui.h"

void DiffMPMLib3D::MaterialPointImGui(MaterialPoint& mp)
{
	ImGui::InputVec3("x", mp.x);
	ImGui::InputVec3("v", mp.v);
	ImGui::InputMat3("F", mp.F);
	ImGui::InputMat3("C", mp.C);
	ImGui::InputMat3("P", mp.P);
	ImGui::InputDouble("m", &mp.m);
	ImGui::InputDouble("vol", &mp.vol);
	ImGui::InputDouble("lam", &mp.lam);
	ImGui::InputDouble("mu", &mp.mu);
	ImGui::InputMat3("dFc", mp.dFc);
}

void DiffMPMLib3D::MaterialPointDisplayImGui(const MaterialPoint& mp)
{
    std::stringstream ss;
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
}
