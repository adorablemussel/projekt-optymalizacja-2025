#pragma once

#include"ode_solver.h"
#include <cmath>
#define _USE_MATH_DEFINES

long long fib_num(int);

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);

matrix ff3T_out(matrix, matrix = NAN, matrix = NAN);
matrix ff3T_in(matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);
matrix df3(double, matrix, matrix = NAN, matrix = NAN);

matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix Hf4T(matrix, matrix = NAN, matrix = NAN);
matrix ff4R(matrix, matrix = NAN, matrix = NAN);
matrix gf4R(matrix, matrix = NAN, matrix = NAN);
double hThetaX(matrix theta, matrix x);

matrix ff5T(matrix, matrix = NAN, matrix = NAN);
matrix ff5R(matrix, matrix = NAN, matrix = NAN);

matrix ff6T(matrix, matrix = NAN, matrix = NAN);
matrix df6(double, matrix, matrix = NAN, matrix = NAN);
matrix ff6R(matrix, matrix = NAN, matrix = NAN);