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

bool ImGui::InputMat3(const char* label, DiffMPMLib3D::Mat3& mat, const char* format, ImGuiInputTextFlags extra_flags)
{
	using namespace DiffMPMLib3D;
	Vec3 row0 = mat.row(0), row1 = mat.row(1), row2 = mat.row(2);
	
	std::string row0label = std::string(label) + "[0,:]";
	bool row0change = InputVec3(row0label.c_str(), row0, format, extra_flags);
	std::string row1label = std::string(label) + "[1,:]";
	bool row1change = InputVec3(row1label.c_str(), row1, format, extra_flags);
	std::string row2label = std::string(label) + "[2,:]";
	bool row2change = InputVec3(row2label.c_str(), row2, format, extra_flags);

	return row0change || row1change || row2change;
}
