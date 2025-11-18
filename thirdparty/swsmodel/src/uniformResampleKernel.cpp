//
// File: uniformResampleKernel.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "uniformResampleKernel.h"
#include "ResampleParser.h"
#include "rt_nonfinite.h"
#include "upfirdnCoreImpl.h"
#include "coder_array.h"
#include <cmath>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : array<double, 2U> &xIn
//                const ResampleParser &opts
// Return Type  : void
//
namespace coder {
namespace b_signal {
namespace internal {
namespace resample {
void uniformResampleAlongFirstDim(array<double, 2U> &xIn,
                                  const ResampleParser &opts)
{
  array<double, 2U> c_xIn;
  array<double, 2U> y;
  array<double, 1U> b_xIn;
  double dv[2];
  if ((!(opts.p == 1.0)) || (!(opts.q == 1.0))) {
    if ((xIn.size(0) == 1) || (xIn.size(1) == 1)) {
      __m128d r;
      double yLen;
      int scalarLB;
      int u1;
      int vectorUB;
      if ((xIn.size(0) == 0) || (xIn.size(1) == 0)) {
        u1 = 0;
      } else {
        scalarLB = xIn.size(0);
        u1 = xIn.size(1);
        if (scalarLB >= u1) {
          u1 = scalarLB;
        }
      }
      yLen = std::ceil(static_cast<double>(u1) * opts.p / opts.q);
      if ((xIn.size(0) == 1) && (xIn.size(1) == 1)) {
        boolean_T xIsRow;
        xIsRow = (xIn.size(0) == 1);
        if (xIsRow) {
          xIn.set_size(xIn.size(0) * xIn.size(1), 1);
        }
        c_xIn.set_size(xIn.size(0), xIn.size(1));
        u1 = xIn.size(0) * xIn.size(1) - 1;
        for (int i{0}; i <= u1; i++) {
          c_xIn[i] = xIn[i];
        }
        upfirdn::upfirdnCoreImpl(
            c_xIn, opts.filterWithPadding, static_cast<double>(xIn.size(0)),
            static_cast<double>(opts.filterWithPadding.size(0)),
            static_cast<double>(xIn.size(1)), opts.p, opts.q, xIn);
        if (xIsRow) {
          u1 = xIn.size(0) * xIn.size(1);
          c_xIn.set_size(1, u1);
          for (int i{0}; i < u1; i++) {
            c_xIn[c_xIn.size(0) * i] = xIn[i];
          }
          xIn.set_size(1, u1);
          for (int i{0}; i < u1; i++) {
            xIn[xIn.size(0) * i] = c_xIn[c_xIn.size(0) * i];
          }
        }
        if ((xIn.size(0) == 0) || (xIn.size(1) == 0)) {
          u1 = 0;
        } else {
          scalarLB = xIn.size(0);
          u1 = xIn.size(1);
          if (scalarLB >= u1) {
            u1 = scalarLB;
          }
        }
        xIn.set_size(u1, 1);
      } else {
        boolean_T xIsRow;
        xIsRow = (xIn.size(0) == 1);
        if (xIsRow) {
          xIn.set_size(xIn.size(0) * xIn.size(1), 1);
        }
        c_xIn.set_size(xIn.size(0), xIn.size(1));
        u1 = xIn.size(0) * xIn.size(1) - 1;
        for (int i{0}; i <= u1; i++) {
          c_xIn[i] = xIn[i];
        }
        upfirdn::upfirdnCoreImpl(
            c_xIn, opts.filterWithPadding, static_cast<double>(xIn.size(0)),
            static_cast<double>(opts.filterWithPadding.size(0)),
            static_cast<double>(xIn.size(1)), opts.p, opts.q, xIn);
        if (xIsRow) {
          u1 = xIn.size(0) * xIn.size(1);
          c_xIn.set_size(1, u1);
          for (int i{0}; i < u1; i++) {
            c_xIn[c_xIn.size(0) * i] = xIn[i];
          }
          xIn.set_size(1, u1);
          for (int i{0}; i < u1; i++) {
            xIn[xIn.size(0) * i] = c_xIn[c_xIn.size(0) * i];
          }
        }
      }
      if (std::isnan(yLen)) {
        y.set_size(1, 1);
        y[0] = rtNaN;
      } else if (yLen < 1.0) {
        y.set_size(y.size(0), 0);
      } else {
        y.set_size(1, static_cast<int>(yLen - 1.0) + 1);
        u1 = static_cast<int>(yLen - 1.0);
        scalarLB = ((static_cast<int>(yLen - 1.0) + 1) / 2) << 1;
        vectorUB = scalarLB - 2;
        for (int i{0}; i <= vectorUB; i += 2) {
          dv[0] = i;
          dv[1] = i + 1;
          r = _mm_loadu_pd(&dv[0]);
          _mm_storeu_pd(&y[i], _mm_add_pd(_mm_set1_pd(1.0), r));
        }
        for (int i{scalarLB}; i <= u1; i++) {
          y[i] = static_cast<double>(i) + 1.0;
        }
      }
      y.set_size(1, y.size(1));
      u1 = y.size(1) - 1;
      scalarLB = (y.size(1) / 2) << 1;
      vectorUB = scalarLB - 2;
      for (int i{0}; i <= vectorUB; i += 2) {
        r = _mm_loadu_pd(&y[i]);
        _mm_storeu_pd(&y[i], _mm_add_pd(_mm_set1_pd(opts.filterDelay), r));
      }
      for (int i{scalarLB}; i <= u1; i++) {
        y[i] = opts.filterDelay + y[i];
      }
      u1 = static_cast<int>(yLen);
      b_xIn.set_size(static_cast<int>(yLen));
      for (int i{0}; i < u1; i++) {
        b_xIn[i] = xIn[static_cast<int>(y[i]) - 1];
      }
      xIn.set_size(static_cast<int>(yLen), 1);
      for (int i{0}; i < u1; i++) {
        xIn[i] = b_xIn[i];
      }
    } else {
      double yLen;
      int scalarLB;
      int u1;
      boolean_T xIsRow;
      yLen = std::ceil(static_cast<double>(xIn.size(0)) * opts.p / opts.q);
      xIsRow = (xIn.size(0) == 1);
      if (xIsRow) {
        xIn.set_size(xIn.size(0) * xIn.size(1), 1);
      }
      c_xIn.set_size(xIn.size(0), xIn.size(1));
      u1 = xIn.size(0) * xIn.size(1) - 1;
      for (int i{0}; i <= u1; i++) {
        c_xIn[i] = xIn[i];
      }
      upfirdn::upfirdnCoreImpl(
          c_xIn, opts.filterWithPadding, static_cast<double>(xIn.size(0)),
          static_cast<double>(opts.filterWithPadding.size(0)),
          static_cast<double>(xIn.size(1)), opts.p, opts.q, xIn);
      if (xIsRow) {
        u1 = xIn.size(0) * xIn.size(1);
        c_xIn.set_size(1, u1);
        for (int i{0}; i < u1; i++) {
          c_xIn[c_xIn.size(0) * i] = xIn[i];
        }
        xIn.set_size(1, u1);
        for (int i{0}; i < u1; i++) {
          xIn[xIn.size(0) * i] = c_xIn[c_xIn.size(0) * i];
        }
      }
      if (std::isnan(yLen)) {
        y.set_size(1, 1);
        y[0] = rtNaN;
      } else if (yLen < 1.0) {
        y.set_size(1, 0);
      } else {
        int vectorUB;
        y.set_size(1, static_cast<int>(yLen - 1.0) + 1);
        u1 = static_cast<int>(yLen - 1.0);
        scalarLB = ((static_cast<int>(yLen - 1.0) + 1) / 2) << 1;
        vectorUB = scalarLB - 2;
        for (int i{0}; i <= vectorUB; i += 2) {
          __m128d r;
          dv[0] = i;
          dv[1] = i + 1;
          r = _mm_loadu_pd(&dv[0]);
          _mm_storeu_pd(&y[i], _mm_add_pd(_mm_set1_pd(1.0), r));
        }
        for (int i{scalarLB}; i <= u1; i++) {
          y[i] = static_cast<double>(i) + 1.0;
        }
      }
      u1 = xIn.size(1);
      scalarLB = y.size(1);
      c_xIn.set_size(y.size(1), u1);
      for (int i{0}; i < u1; i++) {
        for (int i1{0}; i1 < scalarLB; i1++) {
          c_xIn[i1 + c_xIn.size(0) * i] =
              xIn[(static_cast<int>(opts.filterDelay + y[i1]) +
                   xIn.size(0) * i) -
                  1];
        }
      }
      xIn.set_size(y.size(1), xIn.size(1));
      for (int i{0}; i < u1; i++) {
        for (int i1{0}; i1 < scalarLB; i1++) {
          xIn[i1 + xIn.size(0) * i] = c_xIn[i1 + c_xIn.size(0) * i];
        }
      }
    }
  }
}

} // namespace resample
} // namespace internal
} // namespace b_signal
} // namespace coder

//
// File trailer for uniformResampleKernel.cpp
//
// [EOF]
//
