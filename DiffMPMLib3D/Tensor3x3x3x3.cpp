#include "pch.h"
#include "Tensor3x3x3x3.h"


std::string DiffMPMLib3D::ToString(const Tensor3x3x3x3& tensor)
{
    std::streamsize prev_precision = std::cout.precision(16);
    std::stringstream ss;


    auto LeadingSpacesDoubleStr = [](double number, size_t total_width)
    {
        std::string num_str = std::to_string(number);
        std::string ret = "";
        size_t num_chars = num_str.length();
        for (size_t j = 0; j < total_width - num_chars; j++) {
            ret += " ";
        }
        ret += num_str;
        return ret;
    };

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    ss << LeadingSpacesDoubleStr(tensor[a][b](i, j), 20) << ", ";
                }
            }
            ss << std::endl;
        }
    }
    std::cout.precision(prev_precision);
    return ss.str();
}

void DiffMPMLib3D::TensorDiffStats(const Tensor3x3x3x3& A, const Tensor3x3x3x3& B)
{
    // 1. L2 norm of difference
    double l2norm = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double d = A[a][b](i, j) - B[a][b](i, j);
                    l2norm += d * d;
                }
            }
        }
    }
    l2norm = sqrt(l2norm);
    std::cout << "l2 norm = " << l2norm << std::endl;

    // 2. Max difference and index

    double max_diff = 0.0;
    int max_index[4] = { 0,0,0,0 };
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double d = std::fabs(A[a][b](i, j) - B[a][b](i, j));
                    
                    if (d > max_diff) {
                        max_diff = d;
                        max_index[0] = a;
                        max_index[1] = b;
                        max_index[2] = i;
                        max_index[3] = j;
                    }
                }
            }
        }
    }

    std::cout << "max diff = " << max_diff << " at index = " << max_index[0] << ", " << max_index[1] << ", " << max_index[2] << ", " << max_index[3] << std::endl;
}