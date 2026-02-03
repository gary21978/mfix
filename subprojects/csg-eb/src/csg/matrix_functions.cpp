#include "csg_matrix_functions.hpp"

namespace {
double determinant(const matrix::Mat2d &m) {
  return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

double determinant(const matrix::Mat3d &m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) +
         m[0][1] * (m[1][2] * m[2][0] - m[2][2] * m[1][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
}

} // namespace

namespace matrix {

Mat2d inverse(const Mat2d &m) {
  Mat2d res;
  double d = 1 / determinant(m);

  res[0][0] = d * m[1][1];
  res[0][1] = -1.0 * d * m[0][1];
  res[1][0] = -1.0 * d * m[1][0];
  res[1][1] = d * m[0][0];

  return res;
}

Mat3d inverse(const Mat3d &m) {
  Mat3d res;
  double d = 1 / determinant(m);

  res[0][0] = d * (m[1][1] * m[2][2] - m[2][1] * m[1][2]);
  res[0][1] = d * (m[0][2] * m[2][1] - m[0][1] * m[2][2]);
  res[0][2] = d * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
  res[1][0] = d * (m[1][2] * m[2][0] - m[1][0] * m[2][2]);
  res[1][1] = d * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
  res[1][2] = d * (m[1][0] * m[0][2] - m[0][0] * m[1][2]);
  res[2][0] = d * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
  res[2][1] = d * (m[2][0] * m[0][1] - m[0][0] * m[2][1]);
  res[2][2] = d * (m[0][0] * m[1][1] - m[1][0] * m[0][1]);

  return res;
}

} // namespace matrix
