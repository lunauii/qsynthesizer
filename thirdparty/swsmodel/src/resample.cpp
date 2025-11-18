//
// File: resample.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "resample.h"
#include "ResampleParser.h"
#include "firls.h"
#include "kaiser.h"
#include "rt_nonfinite.h"
#include "uniformResampleKernel.h"
#include "upfirdnCoreImpl.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : const array<double, 2U> &varargin_1
//                double varargin_3
//                array<double, 2U> &varargout_1
// Return Type  : void
//
namespace coder {
void resample(const array<double, 2U> &varargin_1, double varargin_3,
              array<double, 2U> &varargout_1)
{
  __m128d r;
  b_signal::internal::resample::ResampleParser opts;
  array<double, 2U> b_varargout_1;
  array<double, 2U> h;
  array<double, 1U> h1;
  array<int, 2U> r2;
  double dv[4];
  double dv1[2];
  double L;
  double absx;
  double nZeroBegin;
  double p1;
  double q1;
  double x;
  double xin;
  int Lx;
  int hi;
  int lastBlockLength;
  int nblocks;
  boolean_T b;
  b = ((varargin_1.size(0) == 1) || (varargin_1.size(1) == 1));
  if (b) {
    opts.dim = 1;
  } else {
    opts.dim = 2;
    if (varargin_1.size(0) != 1) {
      opts.dim = 1;
    }
  }
  opts.inputSize[0] = varargin_1.size(0);
  opts.inputSize[1] = varargin_1.size(1);
  opts.isRowVectorInput = (varargin_1.size(0) == 1);
  if (opts.dim == 1) {
    if (b) {
      if ((varargin_1.size(0) == 0) || (varargin_1.size(1) == 0)) {
        hi = 0;
      } else {
        Lx = varargin_1.size(0);
        hi = varargin_1.size(1);
        if (Lx >= hi) {
          hi = Lx;
        }
      }
      opts.x.set_size(hi, 1);
      for (int k{0}; k < hi; k++) {
        opts.x[k] = varargin_1[k];
      }
    } else {
      opts.x.set_size(varargin_1.size(0), varargin_1.size(1));
      Lx = varargin_1.size(0) * varargin_1.size(1);
      for (int k{0}; k < Lx; k++) {
        opts.x[k] = varargin_1[k];
      }
    }
  } else {
    Lx = varargin_1.size(1);
    hi = varargin_1.size(0);
    opts.x.set_size(varargin_1.size(1), varargin_1.size(0));
    for (int ii{0}; ii < hi; ii++) {
      for (int k{0}; k < Lx; k++) {
        opts.x[k + opts.x.size(0) * ii] =
            varargin_1[ii + varargin_1.size(0) * k];
      }
    }
  }
  x = 6.0 / varargin_3;
  xin = x;
  absx = std::abs(x);
  if (std::isinf(absx) || std::isnan(absx)) {
    absx = rtNaN;
  } else if (absx < 4.4501477170144028E-308) {
    absx = 4.94065645841247E-324;
  } else {
    std::frexp(absx, &lastBlockLength);
    absx = std::ldexp(1.0, lastBlockLength - 53);
  }
  absx = std::fmax(1.0E-12, absx);
  b = std::isnan(x);
  if (std::isinf(x) || b) {
    if (!b) {
      if (x < 0.0) {
        opts.p = -1.0;
      } else {
        opts.p = 1.0;
      }
    } else {
      opts.p = 0.0;
    }
    opts.q = 0.0;
  } else {
    double d;
    double n;
    n = 1.0;
    d = 0.0;
    L = 0.0;
    q1 = 1.0;
    int exitg1;
    do {
      exitg1 = 0;
      p1 = std::round(x);
      if (!std::isinf(x)) {
        x -= p1;
        nZeroBegin = L;
        L = n;
        n = n * p1 + nZeroBegin;
        nZeroBegin = q1;
        q1 = d;
        d = d * p1 + nZeroBegin;
      } else {
        L = n;
        q1 = d;
        n = x;
        d = 0.0;
      }
      if ((x == 0.0) || (std::abs(n / d - xin) <= absx)) {
        exitg1 = 1;
      } else {
        x = 1.0 / x;
      }
    } while (exitg1 == 0);
    if (std::isnan(d)) {
      absx = rtNaN;
    } else if (d < 0.0) {
      absx = -1.0;
    } else {
      absx = (d > 0.0);
    }
    opts.p = n / absx;
    opts.q = std::abs(d);
  }
  absx = std::fmax(opts.p, opts.q);
  L = 20.0 * absx + 1.0;
  dv[0] = 0.0;
  absx = 2.0 * (0.5 / absx);
  dv[1] = absx;
  dv[2] = absx;
  dv[3] = 1.0;
  eFirls(L - 1.0, dv, h, opts.filterWithPadding);
  kaiser(L, opts.filterWithPadding);
  Lx = h.size(1);
  h1.set_size(h.size(1));
  hi = (h.size(1) / 2) << 1;
  lastBlockLength = hi - 2;
  for (int k{0}; k <= lastBlockLength; k += 2) {
    __m128d r1;
    r = _mm_loadu_pd(&h[k]);
    r1 = _mm_loadu_pd(&opts.filterWithPadding[k]);
    _mm_storeu_pd(&h1[k], _mm_mul_pd(r, r1));
  }
  for (int k{hi}; k < Lx; k++) {
    h1[k] = h[k] * opts.filterWithPadding[k];
  }
  if (h1.size(0) == 0) {
    absx = 0.0;
  } else {
    if (h1.size(0) <= 1024) {
      Lx = h1.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      Lx = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(h1.size(0)) >> 10);
      lastBlockLength = h1.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    absx = h1[0];
    for (int k{2}; k <= Lx; k++) {
      absx += h1[k - 1];
    }
    for (int ii{2}; ii <= nblocks; ii++) {
      Lx = (ii - 1) << 10;
      L = h1[Lx];
      if (ii == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (int k{2}; k <= hi; k++) {
        L += h1[(Lx + k) - 1];
      }
      absx += L;
    }
  }
  Lx = h1.size(0);
  opts.filter.set_size(h1.size(0));
  hi = (h1.size(0) / 2) << 1;
  lastBlockLength = hi - 2;
  for (int k{0}; k <= lastBlockLength; k += 2) {
    r = _mm_loadu_pd(&h1[k]);
    _mm_storeu_pd(
        &opts.filter[k],
        _mm_div_pd(_mm_mul_pd(_mm_set1_pd(opts.p), r), _mm_set1_pd(absx)));
  }
  for (int k{hi}; k < Lx; k++) {
    opts.filter[k] = opts.p * h1[k] / absx;
  }
  L = (static_cast<double>(opts.filter.size(0)) - 1.0) / 2.0;
  if (opts.q == 0.0) {
    absx = L;
    if (L == 0.0) {
      absx = 0.0;
    }
  } else if (std::isnan(opts.q)) {
    absx = rtNaN;
  } else if (std::isinf(opts.q)) {
    if (L > 0.0) {
      absx = L;
    } else if (L < 0.0) {
      absx = opts.q;
    } else {
      absx = 0.0;
    }
  } else {
    absx = std::fmod(L, opts.q);
    if (absx == 0.0) {
      absx = opts.q * 0.0;
    } else if (absx < 0.0) {
      absx += opts.q;
    }
  }
  nZeroBegin = std::floor(opts.q - absx);
  opts.filterDelay = std::floor(std::ceil(L + nZeroBegin) / opts.q);
  absx = static_cast<double>(opts.filter.size(0)) + nZeroBegin;
  L = 0.0;
  q1 = opts.q;
  p1 = opts.p;
  Lx = opts.x.size(0);
  while (std::ceil((((static_cast<double>(Lx) - 1.0) * p1 + absx) + L) / q1) -
             opts.filterDelay <
         std::ceil(static_cast<double>(Lx) * p1 / q1)) {
    L++;
  }
  opts.filterWithPadding.set_size(static_cast<int>(absx + L));
  Lx = static_cast<int>(
      (nZeroBegin + static_cast<double>(opts.filter.size(0))) + L);
  if (Lx - 1 >= 0) {
    std::memset(&opts.filterWithPadding[0], 0,
                static_cast<unsigned int>(Lx) * sizeof(double));
  }
  if (opts.filter.size(0) < 1) {
    h.set_size(1, 0);
  } else {
    h.set_size(1, opts.filter.size(0));
    Lx = opts.filter.size(0) - 1;
    hi = (opts.filter.size(0) / 2) << 1;
    lastBlockLength = hi - 2;
    for (int k{0}; k <= lastBlockLength; k += 2) {
      dv1[0] = k;
      dv1[1] = static_cast<double>(k) + 1.0;
      r = _mm_loadu_pd(&dv1[0]);
      _mm_storeu_pd(&h[k], _mm_add_pd(_mm_set1_pd(1.0), r));
    }
    for (int k{hi}; k <= Lx; k++) {
      h[k] = static_cast<double>(k) + 1.0;
    }
  }
  Lx = h.size(1);
  r2.set_size(1, h.size(1));
  if (static_cast<int>(h.size(1) < 1600)) {
    for (int i{0}; i < Lx; i++) {
      r2[i] = static_cast<int>(nZeroBegin + h[i]);
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (int i = 0; i < Lx; i++) {
      r2[i] = static_cast<int>(nZeroBegin + h[i]);
    }
  }
  for (int k{0}; k < Lx; k++) {
    opts.filterWithPadding[r2[k] - 1] = opts.filter[k];
  }
  if (opts.dim == 1) {
    if (opts.isRowVectorInput) {
      varargout_1.set_size(opts.x.size(0), opts.x.size(1));
      Lx = opts.x.size(0) * opts.x.size(1);
      for (int k{0}; k < Lx; k++) {
        varargout_1[k] = opts.x[k];
      }
      b_signal::internal::resample::uniformResampleAlongFirstDim(varargout_1,
                                                                 opts);
      if ((varargout_1.size(0) == 0) || (varargout_1.size(1) == 0)) {
        Lx = 0;
      } else {
        hi = varargout_1.size(0);
        Lx = varargout_1.size(1);
        if (hi >= Lx) {
          Lx = hi;
        }
      }
      varargout_1.set_size(1, Lx);
    } else {
      varargout_1.set_size(opts.x.size(0), opts.x.size(1));
      Lx = opts.x.size(0) * opts.x.size(1);
      for (int k{0}; k < Lx; k++) {
        varargout_1[k] = opts.x[k];
      }
      b_signal::internal::resample::uniformResampleAlongFirstDim(varargout_1,
                                                                 opts);
    }
  } else {
    if (opts.x.size(0) == 1) {
      int i1;
      absx = std::ceil(opts.p / opts.q);
      nblocks = static_cast<int>(absx);
      varargout_1.set_size(static_cast<int>(absx), opts.x.size(1));
      i1 = opts.x.size(1);
      for (int ii{0}; ii < i1; ii++) {
        b_signal::internal::upfirdn::upfirdnCoreImpl(
            opts.x[ii], opts.filterWithPadding,
            static_cast<double>(opts.filterWithPadding.size(0)), opts.p, opts.q,
            h1);
        if (std::isnan(absx)) {
          h.set_size(1, 1);
          h[0] = rtNaN;
        } else if (absx < 1.0) {
          h.set_size(h.size(0), 0);
        } else {
          h.set_size(1, static_cast<int>(absx - 1.0) + 1);
          Lx = static_cast<int>(absx - 1.0);
          lastBlockLength = ((static_cast<int>(absx - 1.0) + 1) / 2) << 1;
          hi = lastBlockLength - 2;
          for (int k{0}; k <= hi; k += 2) {
            dv1[0] = k;
            dv1[1] = k + 1;
            r = _mm_loadu_pd(&dv1[0]);
            _mm_storeu_pd(&h[k], _mm_add_pd(_mm_set1_pd(1.0), r));
          }
          for (int k{lastBlockLength}; k <= Lx; k++) {
            h[k] = static_cast<double>(k) + 1.0;
          }
        }
        h.set_size(1, h.size(1));
        Lx = h.size(1) - 1;
        lastBlockLength = (h.size(1) / 2) << 1;
        hi = lastBlockLength - 2;
        for (int k{0}; k <= hi; k += 2) {
          r = _mm_loadu_pd(&h[k]);
          _mm_storeu_pd(&h[k], _mm_add_pd(_mm_set1_pd(opts.filterDelay), r));
        }
        for (int k{lastBlockLength}; k <= Lx; k++) {
          h[k] = opts.filterDelay + h[k];
        }
        for (int k{0}; k < nblocks; k++) {
          varargout_1[k + varargout_1.size(0) * ii] =
              h1[static_cast<int>(h[k]) - 1];
        }
      }
    } else {
      varargout_1.set_size(opts.x.size(0), opts.x.size(1));
      Lx = opts.x.size(0) * opts.x.size(1);
      for (int k{0}; k < Lx; k++) {
        varargout_1[k] = opts.x[k];
      }
      b_signal::internal::resample::uniformResampleAlongFirstDim(varargout_1,
                                                                 opts);
    }
    if ((opts.x.size(0) == 1) || (opts.x.size(1) == 1)) {
      varargout_1.set_size(static_cast<int>(opts.inputSize[0]),
                           varargout_1.size(0));
    } else {
      Lx = varargout_1.size(1);
      hi = varargout_1.size(0);
      b_varargout_1.set_size(varargout_1.size(1), varargout_1.size(0));
      for (int k{0}; k < hi; k++) {
        for (int ii{0}; ii < Lx; ii++) {
          b_varargout_1[ii + b_varargout_1.size(0) * k] =
              varargout_1[k + varargout_1.size(0) * ii];
        }
      }
      varargout_1.set_size(Lx, hi);
      Lx = b_varargout_1.size(0) * b_varargout_1.size(1);
      for (int k{0}; k < Lx; k++) {
        varargout_1[k] = b_varargout_1[k];
      }
    }
  }
}

} // namespace coder

//
// File trailer for resample.cpp
//
// [EOF]
//
