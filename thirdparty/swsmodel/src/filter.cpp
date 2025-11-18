//
// File: filter.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "filter.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cstring>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : const array<double, 2U> &x
//                array<double, 2U> &y
// Return Type  : void
//
namespace coder {
void filter(const array<double, 2U> &x, array<double, 2U> &y)
{
  static const double dv[2]{1.0, -0.9};
  __m128d r2;
  array<double, 1U> r;
  int b_j;
  int b_k;
  int b_naxpy;
  int b_nx_m_nb;
  int b_scalarLB;
  int loop_ub;
  int nc;
  int nx;
  y.set_size(x.size(0), x.size(1));
  nx = x.size(0);
  nc = x.size(1);
  if (x.size(1) >= 2) {
    int nx_m_nb;
    nx_m_nb = x.size(1);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        r, loop_ub, b_k, b_nx_m_nb, b_scalarLB, b_naxpy, b_j, r2)

    for (int c = 0; c < nx_m_nb; c++) {
      loop_ub = x.size(0);
      r.set_size(x.size(0));
      if (loop_ub - 1 >= 0) {
        std::memset(&r[0], 0,
                    static_cast<unsigned int>(loop_ub) * sizeof(double));
      }
      if (x.size(0) >= 4) {
        for (b_k = 0; b_k < 2; b_k++) {
          b_nx_m_nb = b_k + 1;
          b_scalarLB = ((((loop_ub - b_k) / 2) << 1) + b_k) + 1;
          b_naxpy = b_scalarLB - 2;
          for (b_j = b_nx_m_nb; b_j <= b_naxpy; b_j += 2) {
            r2 = _mm_loadu_pd(&r[b_j - 1]);
            _mm_storeu_pd(
                &r[b_j - 1],
                _mm_add_pd(
                    r2,
                    _mm_mul_pd(
                        _mm_set1_pd(dv[b_k]),
                        _mm_loadu_pd(&x[((b_j - b_k) + x.size(0) * c) - 1]))));
          }
          for (b_j = b_scalarLB; b_j <= loop_ub; b_j++) {
            r[b_j - 1] =
                r[b_j - 1] + dv[b_k] * x[((b_j - b_k) + x.size(0) * c) - 1];
          }
        }
      } else {
        if (x.size(0) > 2) {
          b_nx_m_nb = 0;
        } else {
          b_nx_m_nb = -1;
        }
        for (b_k = 0; b_k <= b_nx_m_nb; b_k++) {
          r2 = _mm_loadu_pd(&r[0]);
          _mm_storeu_pd(&r[0],
                        _mm_add_pd(r2, _mm_mul_pd(_mm_set1_pd(x[x.size(0) * c]),
                                                  _mm_loadu_pd(&dv[0]))));
        }
        b_naxpy = x.size(0) - b_nx_m_nb;
        b_nx_m_nb += 2;
        for (b_k = b_nx_m_nb; b_k <= loop_ub; b_k++) {
          for (b_j = 0; b_j <= b_naxpy - 2; b_j++) {
            b_scalarLB = (b_k + b_j) - 1;
            r[b_scalarLB] =
                r[b_scalarLB] + x[(b_k + x.size(0) * c) - 1] * dv[b_j];
          }
          b_naxpy--;
        }
      }
      b_nx_m_nb = y.size(0);
      for (b_k = 0; b_k < b_nx_m_nb; b_k++) {
        y[b_k + y.size(0) * c] = r[b_k];
      }
    }
  } else {
    int nx_m_nb;
    y.set_size(x.size(0), x.size(1));
    nx_m_nb = x.size(0) * x.size(1);
    for (int k{0}; k < nx_m_nb; k++) {
      y[k] = 0.0;
    }
    for (int b_c{0}; b_c < nc; b_c++) {
      if (nx >= 4) {
        for (int k{0}; k < 2; k++) {
          int naxpy;
          int scalarLB;
          nx_m_nb = k + 1;
          scalarLB = ((((nx - k) / 2) << 1) + k) + 1;
          naxpy = scalarLB - 2;
          for (int j{nx_m_nb}; j <= naxpy; j += 2) {
            __m128d r1;
            r1 = _mm_loadu_pd(&y[j - 1]);
            _mm_storeu_pd(
                &y[j - 1],
                _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(dv[k]),
                                          _mm_loadu_pd(&x[(j - k) - 1]))));
          }
          for (int j{scalarLB}; j <= nx; j++) {
            y[j - 1] = y[j - 1] + dv[k] * x[(j - k) - 1];
          }
        }
      } else {
        int naxpy;
        if (nx > 2) {
          nx_m_nb = 0;
        } else {
          nx_m_nb = -1;
        }
        for (int k{0}; k <= nx_m_nb; k++) {
          __m128d r1;
          r1 = _mm_loadu_pd(&y[0]);
          _mm_storeu_pd(&y[0],
                        _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(x[0]),
                                                  _mm_loadu_pd(&dv[0]))));
        }
        naxpy = nx - nx_m_nb;
        nx_m_nb += 2;
        for (int k{nx_m_nb}; k <= nx; k++) {
          for (int j{0}; j <= naxpy - 2; j++) {
            int scalarLB;
            scalarLB = (k + j) - 1;
            y[scalarLB] = y[scalarLB] + x[k - 1] * dv[j];
          }
          naxpy--;
        }
      }
    }
  }
}

//
// Arguments    : const double b[9]
//                const array<double, 2U> &x
//                array<double, 2U> &y
// Return Type  : void
//
void filter(const double b[9], const array<double, 2U> &x, array<double, 2U> &y)
{
  array<double, 1U> b_b;
  array<double, 1U> b_y1;
  int loop_ub;
  int nx_m_nb;
  loop_ub = x.size(1);
  b_b.set_size(x.size(1));
  for (int k{0}; k < loop_ub; k++) {
    b_b[k] = x[k];
  }
  b_y1.set_size(x.size(1));
  if (loop_ub - 1 >= 0) {
    std::memset(&b_y1[0], 0,
                static_cast<unsigned int>(loop_ub) * sizeof(double));
  }
  if (b_b.size(0) >= 18) {
    for (int k{0}; k < 9; k++) {
      int naxpy;
      int scalarLB;
      nx_m_nb = k + 1;
      scalarLB = ((((loop_ub - k) / 2) << 1) + k) + 1;
      naxpy = scalarLB - 2;
      for (int j{nx_m_nb}; j <= naxpy; j += 2) {
        __m128d r;
        __m128d r1;
        r = _mm_loadu_pd(&b_b[(j - k) - 1]);
        r1 = _mm_loadu_pd(&b_y1[j - 1]);
        _mm_storeu_pd(&b_y1[j - 1],
                      _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(b[k]), r)));
      }
      for (int j{scalarLB}; j <= loop_ub; j++) {
        b_y1[j - 1] = b_y1[j - 1] + b[k] * b_b[(j - k) - 1];
      }
    }
  } else {
    int naxpy;
    if (b_b.size(0) > 9) {
      nx_m_nb = b_b.size(0) - 10;
    } else {
      nx_m_nb = -1;
    }
    for (int k{0}; k <= nx_m_nb; k++) {
      __m128d r;
      __m128d r1;
      r = _mm_loadu_pd(&b_y1[k]);
      r1 = _mm_set1_pd(b_b[k]);
      _mm_storeu_pd(&b_y1[k],
                    _mm_add_pd(r, _mm_mul_pd(r1, _mm_loadu_pd(&b[0]))));
      r = _mm_loadu_pd(&b_y1[k + 2]);
      _mm_storeu_pd(&b_y1[k + 2],
                    _mm_add_pd(r, _mm_mul_pd(r1, _mm_loadu_pd(&b[2]))));
      r = _mm_loadu_pd(&b_y1[k + 4]);
      _mm_storeu_pd(&b_y1[k + 4],
                    _mm_add_pd(r, _mm_mul_pd(r1, _mm_loadu_pd(&b[4]))));
      r = _mm_loadu_pd(&b_y1[k + 6]);
      _mm_storeu_pd(&b_y1[k + 6],
                    _mm_add_pd(r, _mm_mul_pd(r1, _mm_loadu_pd(&b[6]))));
      b_y1[k + 8] = b_y1[k + 8] + b_b[k] * b[8];
    }
    naxpy = b_b.size(0) - nx_m_nb;
    nx_m_nb += 2;
    for (int k{nx_m_nb}; k <= loop_ub; k++) {
      for (int j{0}; j <= naxpy - 2; j++) {
        int scalarLB;
        scalarLB = (k + j) - 1;
        b_y1[scalarLB] = b_y1[scalarLB] + b_b[k - 1] * b[j];
      }
      naxpy--;
    }
  }
  nx_m_nb = b_y1.size(0);
  y.set_size(1, b_y1.size(0));
  for (int k{0}; k < nx_m_nb; k++) {
    y[k] = b_y1[k];
  }
}

} // namespace coder

//
// File trailer for filter.cpp
//
// [EOF]
//
