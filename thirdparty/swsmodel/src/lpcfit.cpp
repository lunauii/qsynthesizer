//
// File: lpcfit.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "lpcfit.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::array<double, 2U> &in2
//                const coder::array<double, 1U> &in3
// Return Type  : void
//
void binary_expand_op(coder::array<double, 2U> &in1,
                      const coder::array<double, 2U> &in2,
                      const coder::array<double, 1U> &in3)
{
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  in1.set_size(1, in1.size(1));
  if (in3.size(0) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(0);
  }
  in1.set_size(in1.size(0), loop_ub);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in3.size(0) != 1);
  if (static_cast<int>(loop_ub < 1600)) {
    for (int i{0}; i < loop_ub; i++) {
      in1[i] = in2[i * stride_0_1] * in3[i * stride_1_1];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (int i = 0; i < loop_ub; i++) {
      in1[i] = in2[i * stride_0_1] * in3[i * stride_1_1];
    }
  }
}

//
// File trailer for lpcfit.cpp
//
// [EOF]
//
