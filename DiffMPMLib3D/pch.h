// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#include "framework.h"

// autodiff include
#include <autodiff/reverse/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <math.h>
#include <iostream>

typedef Eigen::Vector3d Vec3;
typedef Eigen::Matrix3d Mat3;

#endif //PCH_H
