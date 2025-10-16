#pragma once

#include"ode_solver.h"
#include <cmath>
#define _USE_MATH_DEFINES

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
