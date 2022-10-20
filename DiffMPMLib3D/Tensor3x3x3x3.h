#pragma once
#include "pch.h"

typedef std::array<std::array<Mat3, 3>, 3> Tensor3x3x3x3;

std::string ToString(const Tensor3x3x3x3& tensor);

void TensorDiffStats(const Tensor3x3x3x3& A, const Tensor3x3x3x3& B);