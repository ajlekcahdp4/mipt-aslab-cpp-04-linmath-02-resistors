#include "linear_solver.hpp"

#include <range/v3/action.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include "contiguous_matrix.hpp"
#include "equal.hpp"

namespace throttle {

matrix_d nonsingular_solver(matrix_d &&xtnd_matrix) {
  auto cols = xtnd_matrix.cols();
  auto rows = xtnd_matrix.rows();

  xtnd_matrix.convert_to_row_echelon();

  matrix_d res{rows, 1};
  for (matrix_d::size_type i = 0; i < rows; i++) {
    if (throttle::is_roughly_equal(xtnd_matrix[i][i], 0.0)) throw std::runtime_error("Singular matrix provided");
    res[i][0] = xtnd_matrix[i][cols - 1] / xtnd_matrix[i][i];
  }

  return res;
}

matrix_d nonsingular_solver(const matrix_d &coefs, const matrix_d &col) {
  auto rows = coefs.rows();
  auto cols = coefs.cols() + 1;

  if (!coefs.square()) throw std::runtime_error("Non-square matrix provided");
  throttle::linmath::contiguous_matrix<double> xtnd_matrix{rows, cols};

  for (matrix_d::size_type i = 0; i < rows; i++) {
    auto first_row = coefs[i];
    auto second_row = col[i];
    ranges::copy(ranges::views::concat(first_row, second_row), xtnd_matrix.begin() + i * cols);
  }

  return nonsingular_solver(std::move(xtnd_matrix));
}

} // namespace circuits