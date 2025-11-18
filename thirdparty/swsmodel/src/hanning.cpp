//
// File: hanning.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "hanning.h"
#include "rt_nonfinite.h"
#include "swsmodel_rtwutil.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : double varargin_1
//                array<double, 1U> &w
// Return Type  : void
//
namespace coder {
void hanning(double varargin_1, array<double, 1U> &w)
{
  array<double, 2U> r;
  array<double, 1U> b_w;
  double dv[2];
  double n;
  int wt_size;
  boolean_T guard1;
  if (varargin_1 == std::floor(varargin_1)) {
    n = varargin_1;
  } else {
    n = std::round(varargin_1);
  }
  guard1 = false;
  if (n == 0.0) {
    wt_size = 0;
    guard1 = true;
  } else if (n == 1.0) {
    wt_size = 1;
    guard1 = true;
  } else {
    double b;
    if (std::isinf(n)) {
      b = rtNaN;
    } else {
      b = std::fmod(n, 2.0);
    }
    if (b == 0.0) {
      __m128d r1;
      int i1;
      int loop_ub;
      int nx;
      int scalarLB;
      b = n / 2.0;
      if (b < 1.0) {
        r.set_size(1, 0);
      } else {
        r.set_size(1, static_cast<int>(b - 1.0) + 1);
        wt_size = static_cast<int>(b - 1.0);
        scalarLB = ((static_cast<int>(b - 1.0) + 1) / 2) << 1;
        nx = scalarLB - 2;
        for (int i{0}; i <= nx; i += 2) {
          dv[0] = i;
          dv[1] = i + 1;
          r1 = _mm_loadu_pd(&dv[0]);
          _mm_storeu_pd(&r[i], _mm_add_pd(_mm_set1_pd(1.0), r1));
        }
        for (int i{scalarLB}; i <= wt_size; i++) {
          r[i] = static_cast<double>(i) + 1.0;
        }
      }
      wt_size = r.size(1);
      w.set_size(r.size(1));
      scalarLB = (r.size(1) / 2) << 1;
      nx = scalarLB - 2;
      for (int i{0}; i <= nx; i += 2) {
        r1 = _mm_loadu_pd(&r[i]);
        _mm_storeu_pd(
            &w[i], _mm_div_pd(_mm_mul_pd(_mm_set1_pd(6.2831853071795862), r1),
                              _mm_set1_pd(n + 1.0)));
      }
      for (int i{scalarLB}; i < wt_size; i++) {
        w[i] = 6.2831853071795862 * r[i] / (n + 1.0);
      }
      nx = w.size(0);
      if (static_cast<int>(w.size(0) < 1600)) {
        for (int b_k{0}; b_k < nx; b_k++) {
          w[b_k] = std::cos(w[b_k]);
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int b_k = 0; b_k < nx; b_k++) {
          w[b_k] = std::cos(w[b_k]);
        }
      }
      wt_size = (w.size(0) / 2) << 1;
      scalarLB = wt_size - 2;
      for (int i{0}; i <= scalarLB; i += 2) {
        r1 = _mm_loadu_pd(&w[i]);
        _mm_storeu_pd(&w[i], _mm_mul_pd(_mm_set1_pd(0.5),
                                        _mm_sub_pd(_mm_set1_pd(1.0), r1)));
      }
      for (int i{wt_size}; i < nx; i++) {
        w[i] = 0.5 * (1.0 - w[i]);
      }
      if (w.size(0) < 1) {
        nx = 0;
        i1 = 1;
        wt_size = -1;
      } else {
        nx = w.size(0) - 1;
        i1 = -1;
        wt_size = 0;
      }
      scalarLB = div_s32(wt_size - nx, i1);
      loop_ub = (w.size(0) + scalarLB) + 1;
      b_w.set_size(loop_ub);
      wt_size = w.size(0);
      for (int i{0}; i < wt_size; i++) {
        b_w[i] = w[i];
      }
      for (int i{0}; i <= scalarLB; i++) {
        b_w[i + w.size(0)] = w[nx + i1 * i];
      }
      w.set_size(loop_ub);
      for (int i{0}; i < loop_ub; i++) {
        w[i] = b_w[i];
      }
    } else {
      __m128d r1;
      int i1;
      int loop_ub;
      int nx;
      int scalarLB;
      b = (n + 1.0) / 2.0;
      if (b < 1.0) {
        r.set_size(1, 0);
      } else {
        r.set_size(1, static_cast<int>(b - 1.0) + 1);
        wt_size = static_cast<int>(b - 1.0);
        scalarLB = ((static_cast<int>(b - 1.0) + 1) / 2) << 1;
        nx = scalarLB - 2;
        for (int i{0}; i <= nx; i += 2) {
          dv[0] = i;
          dv[1] = i + 1;
          r1 = _mm_loadu_pd(&dv[0]);
          _mm_storeu_pd(&r[i], _mm_add_pd(_mm_set1_pd(1.0), r1));
        }
        for (int i{scalarLB}; i <= wt_size; i++) {
          r[i] = static_cast<double>(i) + 1.0;
        }
      }
      wt_size = r.size(1);
      w.set_size(r.size(1));
      scalarLB = (r.size(1) / 2) << 1;
      nx = scalarLB - 2;
      for (int i{0}; i <= nx; i += 2) {
        r1 = _mm_loadu_pd(&r[i]);
        _mm_storeu_pd(
            &w[i], _mm_div_pd(_mm_mul_pd(_mm_set1_pd(6.2831853071795862), r1),
                              _mm_set1_pd(n + 1.0)));
      }
      for (int i{scalarLB}; i < wt_size; i++) {
        w[i] = 6.2831853071795862 * r[i] / (n + 1.0);
      }
      nx = w.size(0);
      if (static_cast<int>(w.size(0) < 1600)) {
        for (int k{0}; k < nx; k++) {
          w[k] = std::cos(w[k]);
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int k = 0; k < nx; k++) {
          w[k] = std::cos(w[k]);
        }
      }
      wt_size = (w.size(0) / 2) << 1;
      scalarLB = wt_size - 2;
      for (int i{0}; i <= scalarLB; i += 2) {
        r1 = _mm_loadu_pd(&w[i]);
        _mm_storeu_pd(&w[i], _mm_mul_pd(_mm_set1_pd(0.5),
                                        _mm_sub_pd(_mm_set1_pd(1.0), r1)));
      }
      for (int i{wt_size}; i < nx; i++) {
        w[i] = 0.5 * (1.0 - w[i]);
      }
      if (w.size(0) - 1 < 1) {
        nx = 0;
        i1 = 1;
        wt_size = -1;
      } else {
        nx = w.size(0) - 2;
        i1 = -1;
        wt_size = 0;
      }
      scalarLB = div_s32(wt_size - nx, i1);
      loop_ub = (w.size(0) + scalarLB) + 1;
      b_w.set_size(loop_ub);
      wt_size = w.size(0);
      for (int i{0}; i < wt_size; i++) {
        b_w[i] = w[i];
      }
      for (int i{0}; i <= scalarLB; i++) {
        b_w[i + w.size(0)] = w[nx + i1 * i];
      }
      w.set_size(loop_ub);
      for (int i{0}; i < loop_ub; i++) {
        w[i] = b_w[i];
      }
    }
  }
  if (guard1) {
    w.set_size(wt_size);
    for (int i{0}; i < wt_size; i++) {
      w[0] = 1.0;
    }
  }
}

} // namespace coder

//
// File trailer for hanning.cpp
//
// [EOF]
//
