/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tsimmerman.ss@phystech.edu>, <alex.rom23@mail.ru> wrote this file.  As long as you
 * retain this notice you can do whatever you want with this stuff. If we meet
 * some day, and you think this stuff is worth it, you can buy us a beer in
 * return.
 * ----------------------------------------------------------------------------
 */

#pragma once

#include "matrix.hpp"
#include <concepts>

namespace throttle::linmath {

template <std::floating_point T> matrix<T> nonsingular_solver(matrix<T> &&xtnd_matrix) {
  auto cols = xtnd_matrix.cols(), rows = xtnd_matrix.rows();
  xtnd_matrix.convert_to_row_echelon();

  using size_type = typename matrix<T>::size_type;

  matrix<T> res{cols - 1, 1};
  for (size_type i = 0; i < cols - 1; i++) {
    if (is_roughly_equal(xtnd_matrix[i][i], 0.0)) throw std::runtime_error("Singular matrix provided");
    res[i][0] = xtnd_matrix[i][cols - 1] / xtnd_matrix[i][i];
  }

  for (size_type i = cols - 1; i < rows; i++) {
    if (!is_roughly_equal(xtnd_matrix[i][cols - 1], 0.0)) throw std::runtime_error("Singular matrix provided");
  }

  return res;
}

template <std::floating_point T> matrix<T> nonsingular_solver(const matrix<T> &coefs, const matrix<T> &col) {
  auto rows = coefs.rows(), cols = coefs.cols() + 1;

  if (col.cols() != 1) throw std::invalid_argument("A column should be provided");
  if (rows < cols - 1) throw std::runtime_error("System has an infinite number of solutions");

  contiguous_matrix<double> xtnd_matrix{rows, cols};
  using size_type = typename matrix<T>::size_type;

  for (size_type i = 0; i < rows; i++) {
    auto first_row = coefs[i];
    std::copy(first_row.begin(), first_row.end(), xtnd_matrix[i].begin());
    xtnd_matrix[i][cols - 1] = col[i][0];
  }

  return nonsingular_solver(matrix<T>{std::move(xtnd_matrix)});
}

} // namespace throttle::linmath