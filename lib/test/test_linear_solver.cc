#include <gtest/gtest.h>
#include <iostream>

#include "datastructures/vector.hpp"
#include "linmath/linear_solver.hpp"
#include "linmath/matrix.hpp"

namespace linmath = throttle::linmath;

using linear_equation_d = typename throttle::linmath::linear_equation<double>;
using linear_equation_system_d = typename throttle::linmath::linear_equation_system<double>;

TEST(test_linear_solver, test_1) {
  linmath::matrix_d coefs{3, 3, {1, 1, 1, 0, 2, 5, 2, 5, -1}};
  linmath::matrix_d col{3, 1, {6, -4, 27}};
  linmath::matrix_d sol{3, 1, {5, 3, -2}};

  auto res = linmath::nonsingular_solver(coefs, col);
  EXPECT_EQ(res, sol);
}

TEST(test_linear_solver, test_2) {
  linmath::matrix_d coefs{3, 3, {1, 3, -2, 3, 5, 6, 2, 4, 3}};
  linmath::matrix_d col{3, 1, {5, 7, 8}};
  linmath::matrix_d sol{3, 1, {-15, 8, 2}};

  auto res = linmath::nonsingular_solver(coefs, col);
  EXPECT_EQ(res, sol);
}

TEST(test_linear_solver, dependent_system_1) {
  linmath::matrix_d coefs{3, 2, {-1, 2, 2, 3, 1, -2}};
  linmath::matrix_d col{3, 1, {0, 0, 0}};
  linmath::matrix_d sol{2, 1, {0, 0}};
  auto              res = linmath::nonsingular_solver(coefs, col);
  EXPECT_EQ(res, sol);
}

TEST(test_linear_solver, dependent_system_2) {
  linmath::matrix_d coefs{4, 3, {1, 1, 1, 0, -4, -10, 0, 2, 5, 2, 5, -1}};
  linmath::matrix_d col{4, 1, {6, 8, -4, 27}};
  linmath::matrix_d sol{3, 1, {5, 3, -2}};

  auto res = linmath::nonsingular_solver(coefs, col);
  EXPECT_EQ(res, sol);
}

TEST(test_linear_solver, singular_system) {
  linmath::matrix_d coefs{4, 3, {1, 1, 1, 0, 5, -10, 0, 2, 5, 2, 5, -1}};
  linmath::matrix_d col{4, 1, {6, 8, -4, 27}};
  EXPECT_THROW(linmath::nonsingular_solver(coefs, col), std::runtime_error);
}

TEST(test_equation_system, test_1) {
  linear_equation_d        eq1{{1, -1, 7}};
  linear_equation_d        eq2{{3, 2, 16}};
  linear_equation_system_d eqsys{{eq1, eq2}};

  auto res = eqsys.solve();
  EXPECT_TRUE(res.has_value());
  linmath::matrix_d sol{2, 1, {6, -1}};
  EXPECT_EQ(res.value().cols(), sol.cols());
  EXPECT_EQ(res.value(), sol);
}

TEST(test_equation_system, test_2) {
  linear_equation_d        eq1{{1, 1, 1, 6}};
  linear_equation_d        eq2{{0, 2, 5, -4}};
  linear_equation_d        eq3{{2, 5, -1, 27}};
  linear_equation_system_d eqsys{{eq1, eq2, eq3}};
  linmath::matrix_d        sol{3, 1, {5, 3, -2}};
  auto                     res = eqsys.solve();
  EXPECT_EQ(res.value(), sol);
}

TEST(test_equation_system, test_3) {
  linear_equation_d        eq1{{1, 1, 1, 6}};
  linear_equation_d        eq2{{5, 2, 0, -4}};
  linear_equation_d        eq3{{-1, 5, 2, 27}};
  linear_equation_system_d eqsys{{eq1, eq2, eq3}};
  linmath::matrix_d        sol{3, 1, {-2, 3, 5}};
  auto                     res = eqsys.solve();
  EXPECT_EQ(res.value(), sol);
}