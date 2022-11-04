#include "linear_solver.hpp"

#include "contiguous_matrix.hpp"
#include "equal.hpp"

namespace throttle {

matrix_d nonsingular_solver(matrix_d &&xtnd_matrix) {
  auto cols = xtnd_matrix.cols();
  auto rows = xtnd_matrix.rows();
  xtnd_matrix.convert_to_row_echelon();

  matrix_d res{cols - 1, 1};
  for (matrix_d::size_type i = 0; i < cols - 1; i++) {
    if (is_roughly_equal(xtnd_matrix[i][i], 0.0)) throw std::runtime_error("Singular matrix provided");
    res[i][0] = xtnd_matrix[i][cols - 1] / xtnd_matrix[i][i];
  }

  for (matrix_d::size_type i = cols - 1; i < rows; i++) {
    if (!is_roughly_equal(xtnd_matrix[i][cols - 1], 0.0)) throw std::runtime_error("Singular matrix provided");
  }

  return res;
}

matrix_d nonsingular_solver(const matrix_d &coefs, const matrix_d &col) {
  auto rows = coefs.rows();
  auto cols = coefs.cols() + 1;

  if (col.cols() != 1) throw std::invalid_argument("A column should be provided");

  if (rows < cols - 1) throw std::runtime_error("System has an infinite number of solutions");
  throttle::linmath::contiguous_matrix<double> xtnd_matrix{rows, cols};

  for (matrix_d::size_type i = 0; i < rows; i++) {
    auto first_row = coefs[i];

    std::copy(first_row.begin(), first_row.end(), xtnd_matrix[i].begin());
    xtnd_matrix[i][cols - 1] = col[i][0];
  }

  return nonsingular_solver(std::move(xtnd_matrix));
}

} // namespace throttle