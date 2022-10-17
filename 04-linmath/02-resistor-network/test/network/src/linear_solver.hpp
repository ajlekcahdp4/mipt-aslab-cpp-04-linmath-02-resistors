#pragma once

#include "matrix.hpp"

#include <range/v3/action.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace circuits {
using matrix_f = throttle::linmath::matrix<float>;
using matrix_d = throttle::linmath::matrix<double>;

matrix_d nonsingular_solver(matrix_d &&xtnd_matrix) {
  auto cols = xtnd_matrix.cols();
  auto rows = xtnd_matrix.rows();
  xtnd_matrix.convert_to_row_echelon();
  matrix_d res{rows, 1};
  for (unsigned i = 0; i < rows; i++)
    res[i][0] = xtnd_matrix[i][cols - 1];
  return res;
}

matrix_d linear_solver(const matrix_d &coefs, const matrix_d &col) {
  auto rows = coefs.rows();
  auto cols = coefs.cols() + 1;

  matrix_d xtnd_matrix{rows, cols};
  for (unsigned i = 0; i < rows; i++) {
    auto concated = ranges::views::concat(coefs[i], col[i]);
    ranges::copy(concated, xtnd_matrix.begin());
  }
  return nonsingular_solver(std::move(xtnd_matrix));
}

} // namespace circuits