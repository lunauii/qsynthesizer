//
// File: abs.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "abs.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const creal_T x_data[]
//                const int x_size[2]
//                double y_data[]
//                int y_size[2]
// Return Type  : void
//
namespace coder {
void b_abs(const creal_T x_data[], const int x_size[2], double y_data[],
           int y_size[2])
{
  int i;
  i = x_size[0] * x_size[1];
  y_size[0] = static_cast<signed char>(x_size[0]);
  y_size[1] = static_cast<signed char>(x_size[1]);
  for (int k{0}; k < i; k++) {
    double a;
    double b;
    a = std::abs(x_data[k].re);
    b = std::abs(x_data[k].im);
    if (a < b) {
      a /= b;
      y_data[k] = b * std::sqrt(a * a + 1.0);
    } else if (a > b) {
      b /= a;
      y_data[k] = a * std::sqrt(b * b + 1.0);
    } else if (std::isnan(b)) {
      y_data[k] = rtNaN;
    } else {
      y_data[k] = a * 1.4142135623730951;
    }
  }
}

} // namespace coder

//
// File trailer for abs.cpp
//
// [EOF]
//
