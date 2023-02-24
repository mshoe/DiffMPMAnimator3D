#pragma once
#include "pch.h"
#include "imgui.h"
#include "misc/cpp/imgui_stdlib.h"

namespace ImGui
{
	

	bool InputDouble3(const char* label, double v[3], const char* format = "%.15f", ImGuiInputTextFlags extra_flags = 0);
	bool InputVec3(const char* label, DiffMPMLib3D::Vec3& vec, const char* format = "%.15f", ImGuiInputTextFlags extra_flags = 0);

	bool InputMat3(const char* label, DiffMPMLib3D::Mat3& mat, const char* format = "%.15f", ImGuiInputTextFlags extra_flags = 0);
}