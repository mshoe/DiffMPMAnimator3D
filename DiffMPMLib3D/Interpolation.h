#pragma once
#include "pch.h"

inline double CubicBSpline(double x) {
	x = abs(x);
	if (0.0 <= x && x < 1.0) {
		return 0.5 * x * x * x - x * x + 2.0 / 3.0;
	}
	else if (1.0 <= x && x < 2.0) {
		return (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
	}
	else {
		return 0.0;
	}
}

inline double CubicBSplineSlope(double x)
{
	double absx = abs(x);
	if (0.0 <= absx && absx < 1.0) {
		return 1.5 * x * absx - 2.0 * x;
	}
	else if (1.0 <= absx && absx < 2.0) {
		return -x * absx / 2.0 + 2.0 * x - 2.0 * x / absx;
	}
	else {
		return 0.0;
	}
}