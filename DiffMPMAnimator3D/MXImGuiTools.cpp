#include "MXImGuiTools.h"


bool ImGui::InputDouble3(const char* label, double v[3], const char* format, ImGuiInputTextFlags extra_flags)
{
    return InputScalarN(label, ImGuiDataType_Double, v, 3, NULL, NULL, format, extra_flags);
}

bool ImGui::InputVec3(const char* label, DiffMPMLib3D::Vec3& vec, const char* format, ImGuiInputTextFlags extra_flags)
{
	double vec_array[3]{ vec.x(), vec.y(), vec.z()};
	bool ret = ImGui::InputDouble3(label, vec_array, format, extra_flags);
	vec = DiffMPMLib3D::Vec3(vec_array[0], vec_array[1], vec_array[2]);
	return ret;
}
