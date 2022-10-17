#pragma once

#include "matrix.hpp"

namespace throttle {
using matrix_d = throttle::linmath::matrix<double>;

matrix_d nonsingular_solver(matrix_d &&xtnd_matrix);
matrix_d nonsingular_solver(const matrix_d &coefs, const matrix_d &col);

} // namespace circuits