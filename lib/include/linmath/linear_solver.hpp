/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICEn_eqsSE" (Revision 42):
 * <tsimmerman.ss@phystech.edu>, <alex.rom23@mail.ru> wrote this file.  As long as you
 * retain this notice you can do whatever you want with this stuff. If we meet
 * some day, and you think this stuff is worth it, you can buy us a beer in
 * return.
 * ----------------------------------------------------------------------------
 */

#pragma once

#include "linmath/matrix.hpp"
#include "datastructures/vector.hpp"
#include <concepts>
#include <optional>

namespace throttle::linmath {

template <std::floating_point T> struct linear_equation final {
  using value_type = T;
  using size_type = std::size_t;

  size_type                      n_vars{};
  containers::vector<value_type> m_coefs;

  // if coefs includes free coef
  linear_equation(const containers::vector<value_type> &coefs) : m_coefs{coefs}, n_vars{coefs.size() - 1} {}

  // if free coef provided separately
  linear_equation(const containers::vector<value_type> &coefs, const value_type &free)
      : n_vars{coefs.size()}, m_coefs{n_vars + 1} {
    for (auto i = 0; i < m_coefs.size(); ++i)
      m_coefs[0][i] = coefs[0][i];
    m_coefs[0][m_coefs.size()] = free;
  }

  template <std::input_iterator it>
  linear_equation(it start, it finish) : n_vars{std::distance(start, finish) - 1}, m_coefs{n_vars + 1} {
    std::copy(start, finish, m_coefs.begin());
  }

  linear_equation(std::initializer_list<value_type> list) : linear_equation{list.begin(), list.end()} {}
};

template <std::floating_point T> struct linear_equation_system final {
  using value_type = T;
  using size_type = std::size_t;
  using equation_t = linear_equation<value_type>;

  size_type                                n_eqs;
  throttle::containers::vector<equation_t> equations;

  linear_equation_system(const containers::vector<equation_t> equations) : equations{equations} {}

  matrix<value_type> get_xtnd_matrix() const {
    auto rows = n_eqs;
    auto cols = equations[0].n_vars + 1;

    matrix<value_type> xtnd_matrix{rows, cols};
    for (unsigned i = 0; i < rows; ++i)
      for (unsigned j = 0; j < cols; ++j)
        xtnd_matrix[i][j] = equations[i].m_coefs[j];
    return xtnd_matrix;
  }

  template <std::input_iterator it>
  linear_equation_system(it start, it finish) : n_eqs{std::distance(start, finish)}, equations{n_eqs} {
    std::copy(start, finish, equations.begin());
  }

  linear_equation_system(std::initializer_list<equation_t> list) : linear_equation_system{list.begin(), list.end()} {}

  std::optional<matrix<value_type>> solve() const {
    auto xtnd_matrix = get_xtnd_matrix();
    auto cols = xtnd_matrix.cols();
    auto rows = xtnd_matrix.rows();

    xtnd_matrix.convert_to_row_echelon();

    matrix<T> res{cols - 1, 1};
    for (size_type i = 0; i < cols - 1; i++) {
      if (is_roughly_equal(xtnd_matrix[i][i], 0.0)) return std::nullopt;
      res[i][0] = xtnd_matrix[i][cols - 1] / xtnd_matrix[i][i];
    }

    for (size_type i = cols - 1; i < rows; i++) {
      if (!is_roughly_equal(xtnd_matrix[i][cols - 1], 0.0)) return std::nullopt;
    }

    return res;
  }
};

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