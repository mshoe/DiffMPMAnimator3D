// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#include "framework.h"

#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <math.h>
#include <iostream>

typedef Eigen::Vector3d Vec3;
typedef Eigen::Matrix3d Mat3;

inline double InnerProduct(const Mat3& A, const Mat3& B)
{
	double ret = 0.0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			ret += A(i, j) * B(i, j);
		}
	}
	return ret;
}

inline Mat3 CofactorMatrix(const Mat3& A)
{
	// https://en.wikipedia.org/wiki/Adjugate_matrix
	return A.adjoint().transpose();
}

#endif //PCH_H
