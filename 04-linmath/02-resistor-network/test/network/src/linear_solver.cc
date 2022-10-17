#include "linear_solver.hpp"

namespace circuits {
inline matrix_d nonsingular_solver(matrix_d &&xtnd_matrix) {
  auto cols = xtnd_matrix.cols();
  auto rows = xtnd_matrix.rows();
  xtnd_matrix.convert_to_row_echelon();

  matrix_d res{rows, 1};
  for (unsigned i = 0; i < rows; i++) {
    if (is_roughly_equal(xtnd_matrix[i][i], 0)) {
      if (!is_roughly_equal(xtnd_matrix[i][cols - 1], 0)) throw std::runtime_error("singular matrix provided");
      return 0;
    } else
      res[i][0] = xtnd_matrix[i][cols - 1] / xtnd_matrix[i][i];
  }
  return res;
}

inline matrix_d linear_solver(const matrix_d &coefs, const matrix_d &col) {
  auto rows = coefs.rows();
  auto cols = coefs.cols() + 1;
  matrix_d xtnd_matrix{rows, cols};
  for (unsigned i = 0; i < rows; i++) {
    auto first_row = coefs[i];
    auto second_row = col[i];
    auto concated = ranges::views::concat(first_row, second_row);
    ranges::copy(concated, xtnd_matrix.begin());
  }
  return nonsingular_solver(std::move(xtnd_matrix));
}
} // namespace circuits