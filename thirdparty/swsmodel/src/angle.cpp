//
// File: angle.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "angle.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
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
void angle(const creal_T x_data[], const int x_size[2], double y_data[],
           int y_size[2])
{
  int i;
  i = x_size[0] * x_size[1];
  y_size[0] = static_cast<signed char>(x_size[0]);
  y_size[1] = static_cast<signed char>(x_size[1]);
  for (int k{0}; k < i; k++) {
    double x;
    double y;
    y = x_data[k].im;
    x = x_data[k].re;
    if (std::isnan(y) || std::isnan(x)) {
      y = rtNaN;
    } else if (std::isinf(y) && std::isinf(x)) {
      int i1;
      int i2;
      if (y > 0.0) {
        i1 = 1;
      } else {
        i1 = -1;
      }
      if (x > 0.0) {
        i2 = 1;
      } else {
        i2 = -1;
      }
      y = std::atan2(static_cast<double>(i1), static_cast<double>(i2));
    } else if (x == 0.0) {
      if (y > 0.0) {
        y = RT_PI / 2.0;
      } else if (y < 0.0) {
        y = -(RT_PI / 2.0);
      } else {
        y = 0.0;
      }
    } else {
      y = std::atan2(y, x);
    }
    y_data[k] = y;
  }
}

} // namespace coder

//
// File trailer for angle.cpp
//
// [EOF]
//
