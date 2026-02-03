#ifndef MATRIX_FUNCTIONS_H_
#define MATRIX_FUNCTIONS_H_

#include <array>

namespace matrix {
using Mat2d = std::array<std::array<double, 2>, 2>;
using Mat3d = std::array<std::array<double, 3>, 3>;

Mat2d inverse(const Mat2d &m);
Mat3d inverse(const Mat3d &m);
} // namespace matrix

#endif
