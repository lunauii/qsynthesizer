//
// File: xcorr.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "xcorr.h"
#include "fft.h"
#include "ifft.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "coder_array.h"
#include "omp.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : const array<double, 2U> &x
//                array<double, 2U> &c
// Return Type  : void
//
namespace coder {
void xcorr(const array<double, 2U> &x, array<double, 2U> &c)
{
  array<creal_T, 1U> X;
  array<double, 1U> b_c1;
  array<double, 1U> b_x;
  array<double, 1U> c1;
  array<int, 2U> y;
  double a;
  double b;
  double f_tmp;
  double m2;
  double s;
  double varargin_1;
  int ceilLog2;
  int mxl;
  int nx;
  boolean_T doAutoCorrTd;
  mxl = static_cast<int>(std::fmin(static_cast<double>(x.size(1)) - 1.0,
                                   static_cast<double>(x.size(1)) - 1.0));
  f_tmp = 2.0 * static_cast<double>(x.size(1)) - 1.0;
  s = std::frexp(std::abs(f_tmp), &ceilLog2);
  if (s == 0.5) {
    ceilLog2--;
  }
  nx = static_cast<int>(std::abs(static_cast<double>(ceilLog2)));
  if (nx == 0) {
    m2 = 1.0;
  } else if (nx == 1) {
    if (ceilLog2 > 0) {
      m2 = 2.0;
    } else {
      m2 = 0.5;
    }
  } else if (ceilLog2 == 2) {
    m2 = 4.0;
  } else {
    m2 = std::pow(2.0, static_cast<double>(ceilLog2));
  }
  nx = x.size(1);
  doAutoCorrTd = true;
  for (int k{0}; k < nx; k++) {
    if (doAutoCorrTd) {
      s = x[k];
      if (std::isinf(s) || std::isnan(s)) {
        doAutoCorrTd = false;
      }
    } else {
      doAutoCorrTd = false;
    }
  }
  doAutoCorrTd = ((!doAutoCorrTd) ||
                  (f_tmp + static_cast<double>(mxl) *
                               ((f_tmp - static_cast<double>(mxl)) - 1.0) <
                   m2 * (15.0 * static_cast<double>(ceilLog2) + 6.0)));
  if (doAutoCorrTd) {
    int b_mxl;
    b_mxl = static_cast<int>(std::fmin(static_cast<double>(mxl),
                                       static_cast<double>(x.size(1)) - 1.0));
    nx = static_cast<int>(2.0 * static_cast<double>(b_mxl) + 1.0);
    c1.set_size(nx);
    if (nx - 1 >= 0) {
      std::memset(&c1[0], 0, static_cast<unsigned int>(nx) * sizeof(double));
    }
    for (int k{0}; k <= b_mxl; k++) {
      s = 0.0;
      ceilLog2 = x.size(1) - k;
      for (int b_i{0}; b_i < ceilLog2; b_i++) {
        s += x[b_i] * x[k + b_i];
      }
      c1[b_mxl - k] = s;
      c1[b_mxl + k] = s;
    }
  } else {
    int b_mxl;
    b_x = x.reshape(x.size(1));
    fft(b_x, m2, X);
    ceilLog2 = X.size(0);
    b_c1.set_size(X.size(0));
    nx = (X.size(0) < 1600);
    if (nx) {
      for (int b_k{0}; b_k < ceilLog2; b_k++) {
        a = std::abs(X[b_k].re);
        b = std::abs(X[b_k].im);
        if (a < b) {
          a /= b;
          b_c1[b_k] = b * std::sqrt(a * a + 1.0);
        } else if (a > b) {
          b /= a;
          b_c1[b_k] = a * std::sqrt(b * b + 1.0);
        } else if (std::isnan(b)) {
          b_c1[b_k] = rtNaN;
        } else {
          b_c1[b_k] = a * 1.4142135623730951;
        }
      }
      for (int i{0}; i < ceilLog2; i++) {
        s = b_c1[i];
        b_c1[i] = s * s;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b, a)

      for (int b_k = 0; b_k < ceilLog2; b_k++) {
        a = std::abs(X[b_k].re);
        b = std::abs(X[b_k].im);
        if (a < b) {
          a /= b;
          b_c1[b_k] = b * std::sqrt(a * a + 1.0);
        } else if (a > b) {
          b /= a;
          b_c1[b_k] = a * std::sqrt(b * b + 1.0);
        } else if (std::isnan(b)) {
          b_c1[b_k] = rtNaN;
        } else {
          b_c1[b_k] = a * 1.4142135623730951;
        }
      }
#pragma omp parallel for num_threads(omp_get_max_threads()) private(varargin_1)

      for (int i = 0; i < ceilLog2; i++) {
        varargin_1 = b_c1[i];
        b_c1[i] = varargin_1 * varargin_1;
      }
    }
    ifft(b_c1, X);
    nx = X.size(0);
    b_c1.set_size(X.size(0));
    for (int k{0}; k < nx; k++) {
      b_c1[k] = X[k].re;
    }
    if (mxl < 0) {
      b_mxl = -1;
    } else {
      b_mxl = mxl;
    }
    if (mxl < 1) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, mxl);
      nx = (mxl / 4) << 2;
      ceilLog2 = nx - 4;
      for (int k{0}; k <= ceilLog2; k += 4) {
        _mm_storeu_si128(
            (__m128i *)&y[k],
            _mm_add_epi32(
                _mm_set1_epi32(1),
                _mm_add_epi32(_mm_set1_epi32(k),
                              _mm_loadu_si128((const __m128i *)&iv[0]))));
      }
      for (int k{nx}; k < mxl; k++) {
        y[k] = k + 1;
      }
    }
    s = m2 - static_cast<double>(mxl);
    c1.set_size((y.size(1) + b_mxl) + 1);
    nx = y.size(1);
    for (int k{0}; k < nx; k++) {
      c1[k] = b_c1[static_cast<int>(s + static_cast<double>(y[k])) - 1];
    }
    for (int k{0}; k <= b_mxl; k++) {
      c1[k + y.size(1)] = b_c1[k];
    }
  }
  ceilLog2 =
      static_cast<int>(2.0 * (static_cast<double>(x.size(1)) - 1.0) + 1.0);
  b_c1.set_size(ceilLog2);
  if (ceilLog2 - 1 >= 0) {
    std::memset(&b_c1[0], 0,
                static_cast<unsigned int>(ceilLog2) * sizeof(double));
  }
  nx = (x.size(1) - 1) << 1;
  if (nx >= 0) {
    std::copy(&c1[0], &c1[nx + 1], &b_c1[0]);
  }
  c.set_size(1, ceilLog2);
  for (int k{0}; k < ceilLog2; k++) {
    c[k] = b_c1[k];
  }
}

} // namespace coder

//
// File trailer for xcorr.cpp
//
// [EOF]
//
