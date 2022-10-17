#pragma once

#include "equal.hpp"
#include "matrix.hpp"

#include <range/v3/action.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace circuits {
using matrix_f = throttle::linmath::matrix<float>;
using matrix_d = throttle::linmath::matrix<double>;

inline matrix_d nonsingular_solver(matrix_d &&xtnd_matrix);

inline matrix_d linear_solver(const matrix_d &coefs, const matrix_d &col);

} // namespace circuits