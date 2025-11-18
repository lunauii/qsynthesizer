//
// File: padarray.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "padarray.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
//
// Arguments    : const array<double, 2U> &varargin_1
//                double varargin_2
//                array<double, 2U> &b
// Return Type  : void
//
namespace coder {
void padarray(const array<double, 2U> &varargin_1, double varargin_2,
              array<double, 2U> &b)
{
  if ((varargin_1.size(0) == 0) || (varargin_1.size(1) == 0)) {
    double padSize_idx_0;
    int loop_ub;
    padSize_idx_0 = static_cast<double>(varargin_1.size(0)) + 2.0 * varargin_2;
    b.set_size(static_cast<int>(padSize_idx_0), varargin_1.size(1));
    loop_ub = static_cast<int>(padSize_idx_0) * varargin_1.size(1);
    for (int j{0}; j < loop_ub; j++) {
      b[j] = 0.0;
    }
  } else {
    int i;
    int i1;
    int loop_ub;
    i = static_cast<int>(static_cast<double>(varargin_1.size(0)) +
                         2.0 * varargin_2);
    loop_ub = varargin_1.size(1);
    b.set_size(i, varargin_1.size(1));
    i1 = varargin_1.size(1) + 1;
    for (int j{i1}; j <= loop_ub; j++) {
      for (int b_i{0}; b_i < i; b_i++) {
        b[b_i + b.size(0) * (j - 1)] = 0.0;
      }
    }
    i1 = varargin_1.size(1);
    loop_ub = static_cast<int>(varargin_2);
    for (int j{0}; j < i1; j++) {
      for (int b_i{0}; b_i < loop_ub; b_i++) {
        b[b_i + b.size(0) * j] = 0.0;
      }
    }
    loop_ub = (static_cast<int>(varargin_2) + varargin_1.size(0)) + 1;
    for (int j{0}; j < i1; j++) {
      for (int b_i{loop_ub}; b_i <= i; b_i++) {
        b[(b_i + b.size(0) * j) - 1] = 0.0;
      }
    }
    loop_ub = varargin_1.size(0);
    for (int j{0}; j < i1; j++) {
      for (int b_i{0}; b_i < loop_ub; b_i++) {
        b[(b_i + static_cast<int>(varargin_2)) + b.size(0) * j] =
            varargin_1[b_i + varargin_1.size(0) * j];
      }
    }
  }
}

} // namespace coder

//
// File trailer for padarray.cpp
//
// [EOF]
//
