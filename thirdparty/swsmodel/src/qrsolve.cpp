//
// File: qrsolve.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "qrsolve.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "xnrm2.h"
#include "xzlarfg.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : const array<double, 2U> &A
//                const array<double, 1U> &B
//                array<double, 1U> &Y
// Return Type  : void
//
namespace coder {
namespace internal {
void qrsolve(const array<double, 2U> &A, const array<double, 1U> &B,
             array<double, 1U> &Y)
{
  array<double, 2U> b_A;
  array<double, 1U> tau;
  array<double, 1U> vn1;
  array<double, 1U> vn2;
  array<double, 1U> work;
  array<int, 2U> jpvt;
  double d;
  double tol;
  int b_m;
  int m;
  int ma;
  int maxmn;
  int minmana;
  int minmn;
  int mn;
  int n;
  int rankA;
  m = A.size(0);
  b_m = A.size(1);
  b_A.set_size(A.size(0), A.size(1));
  minmn = A.size(0) * A.size(1);
  for (int j{0}; j < minmn; j++) {
    b_A[j] = A[j];
  }
  ma = A.size(0);
  n = A.size(1) - 1;
  minmn = A.size(0);
  minmana = A.size(1);
  if (minmn <= minmana) {
    minmana = minmn;
  }
  tau.set_size(minmana);
  if (minmana - 1 >= 0) {
    std::memset(&tau[0], 0,
                static_cast<unsigned int>(minmana) * sizeof(double));
  }
  if ((A.size(0) == 0) || (A.size(1) == 0)) {
    jpvt.set_size(1, A.size(1));
    if (b_m - 1 >= 0) {
      std::memset(&jpvt[0], 0, static_cast<unsigned int>(b_m) * sizeof(int));
    }
    minmn = (A.size(1) / 4) << 2;
    maxmn = minmn - 4;
    for (int j{0}; j <= maxmn; j += 4) {
      _mm_storeu_si128(
          (__m128i *)&jpvt[j],
          _mm_add_epi32(_mm_add_epi32(_mm_set1_epi32(j),
                                      _mm_loadu_si128((const __m128i *)&iv[0])),
                        _mm_set1_epi32(1)));
    }
    for (int j{minmn}; j <= n; j++) {
      jpvt[j] = j + 1;
    }
  } else {
    jpvt.set_size(1, A.size(1));
    std::memset(&jpvt[0], 0, static_cast<unsigned int>(b_m) * sizeof(int));
    minmn = (A.size(1) / 4) << 2;
    maxmn = minmn - 4;
    for (int j{0}; j <= maxmn; j += 4) {
      _mm_storeu_si128(
          (__m128i *)&jpvt[j],
          _mm_add_epi32(_mm_add_epi32(_mm_set1_epi32(j),
                                      _mm_loadu_si128((const __m128i *)&iv[0])),
                        _mm_set1_epi32(1)));
    }
    for (int j{minmn}; j <= n; j++) {
      jpvt[j] = j + 1;
    }
    work.set_size(A.size(1));
    vn1.set_size(A.size(1));
    vn2.set_size(A.size(1));
    std::memset(&work[0], 0, static_cast<unsigned int>(b_m) * sizeof(double));
    std::memset(&vn1[0], 0, static_cast<unsigned int>(b_m) * sizeof(double));
    std::memset(&vn2[0], 0, static_cast<unsigned int>(b_m) * sizeof(double));
    if (static_cast<int>(A.size(1) < 1600)) {
      for (int b_j{0}; b_j <= n; b_j++) {
        d = blas::xnrm2(m, A, b_j * m + 1);
        vn1[b_j] = d;
        vn2[b_j] = d;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(d)

      for (int b_j = 0; b_j <= n; b_j++) {
        d = blas::xnrm2(ma, A, b_j * ma + 1);
        vn1[b_j] = d;
        vn2[b_j] = d;
      }
    }
    for (int i{0}; i < minmana; i++) {
      double atmp;
      double temp1;
      int ii;
      int ip1;
      int mmi;
      ip1 = i + 2;
      mn = i * ma;
      ii = mn + i;
      rankA = n - i;
      mmi = (ma - i) - 1;
      minmn = rankA + 1;
      if (rankA < 0) {
        maxmn = -1;
      } else {
        maxmn = 0;
        if (rankA > 0) {
          tol = std::abs(vn1[i]);
          for (int j{2}; j <= minmn; j++) {
            temp1 = std::abs(vn1[(i + j) - 1]);
            if (temp1 > tol) {
              maxmn = j - 1;
              tol = temp1;
            }
          }
        }
      }
      b_m = i + maxmn;
      if (b_m + 1 != i + 1) {
        minmn = b_m * ma;
        for (int j{0}; j < ma; j++) {
          maxmn = minmn + j;
          tol = b_A[maxmn];
          m = mn + j;
          b_A[maxmn] = b_A[m];
          b_A[m] = tol;
        }
        minmn = jpvt[b_m];
        jpvt[b_m] = jpvt[i];
        jpvt[i] = minmn;
        vn1[b_m] = vn1[i];
        vn2[b_m] = vn2[i];
      }
      if (i + 1 < ma) {
        tol = b_A[ii];
        temp1 = reflapack::xzlarfg(mmi + 1, tol, b_A, ii + 2);
        tau[i] = temp1;
        b_A[ii] = tol;
      } else {
        temp1 = 0.0;
        tau[i] = 0.0;
      }
      if (i < n) {
        int lastv;
        atmp = b_A[ii];
        b_A[ii] = 1.0;
        mn = (ii + ma) + 1;
        if (temp1 != 0.0) {
          boolean_T exitg2;
          lastv = mmi;
          minmn = ii + mmi;
          while ((lastv + 1 > 0) && (b_A[minmn] == 0.0)) {
            lastv--;
            minmn--;
          }
          b_m = rankA - 1;
          exitg2 = false;
          while ((!exitg2) && (b_m + 1 > 0)) {
            int exitg1;
            minmn = mn + b_m * ma;
            maxmn = minmn;
            do {
              exitg1 = 0;
              if (maxmn <= minmn + lastv) {
                if (b_A[maxmn - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  maxmn++;
                }
              } else {
                b_m--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);
            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = -1;
          b_m = -1;
        }
        if (lastv + 1 > 0) {
          if (b_m + 1 != 0) {
            if (b_m >= 0) {
              std::memset(&work[0], 0,
                          static_cast<unsigned int>(b_m + 1) * sizeof(double));
            }
            minmn = 0;
            maxmn = mn + ma * b_m;
            for (int j{mn}; ma < 0 ? j >= maxmn : j <= maxmn; j += ma) {
              tol = 0.0;
              m = j + lastv;
              for (int c_j{j}; c_j <= m; c_j++) {
                tol += b_A[c_j - 1] * b_A[(ii + c_j) - j];
              }
              work[minmn] = work[minmn] + tol;
              minmn++;
            }
          }
          if (!(-tau[i] == 0.0)) {
            for (int j{0}; j <= b_m; j++) {
              if (work[j] != 0.0) {
                tol = work[j] * -tau[i];
                minmn = lastv + mn;
                for (int c_j{mn}; c_j <= minmn; c_j++) {
                  b_A[c_j - 1] = b_A[c_j - 1] + b_A[(ii + c_j) - mn] * tol;
                }
              }
              mn += ma;
            }
          }
        }
        b_A[ii] = atmp;
      }
      for (int j{ip1}; j <= n + 1; j++) {
        minmn = i + (j - 1) * ma;
        tol = vn1[j - 1];
        if (tol != 0.0) {
          temp1 = std::abs(b_A[minmn]) / tol;
          temp1 = 1.0 - temp1 * temp1;
          if (temp1 < 0.0) {
            temp1 = 0.0;
          }
          atmp = tol / vn2[j - 1];
          atmp = temp1 * (atmp * atmp);
          if (atmp <= 1.4901161193847656E-8) {
            if (i + 1 < ma) {
              tol = blas::xnrm2(mmi, b_A, minmn + 2);
              vn1[j - 1] = tol;
              vn2[j - 1] = tol;
            } else {
              vn1[j - 1] = 0.0;
              vn2[j - 1] = 0.0;
            }
          } else {
            vn1[j - 1] = tol * std::sqrt(temp1);
          }
        }
      }
    }
  }
  rankA = 0;
  if (b_A.size(0) < b_A.size(1)) {
    minmn = b_A.size(0);
    maxmn = b_A.size(1);
  } else {
    minmn = b_A.size(1);
    maxmn = b_A.size(0);
  }
  if (minmn > 0) {
    tol = std::fmin(1.4901161193847656E-8,
                    2.2204460492503131E-15 * static_cast<double>(maxmn)) *
          std::abs(b_A[0]);
    while ((rankA < minmn) &&
           (!(std::abs(b_A[rankA + b_A.size(0) * rankA]) <= tol))) {
      rankA++;
    }
  }
  minmn = B.size(0);
  work.set_size(B.size(0));
  for (int j{0}; j < minmn; j++) {
    work[j] = B[j];
  }
  minmn = b_A.size(1);
  Y.set_size(b_A.size(1));
  for (int j{0}; j < minmn; j++) {
    Y[j] = 0.0;
  }
  minmn = b_A.size(0);
  mn = b_A.size(1);
  if (minmn <= mn) {
    mn = minmn;
  }
  for (int c_j{0}; c_j < mn; c_j++) {
    b_m = b_A.size(0);
    if (tau[c_j] != 0.0) {
      tol = work[c_j];
      minmn = c_j + 2;
      for (int j{minmn}; j <= b_m; j++) {
        tol += b_A[(j + b_A.size(0) * c_j) - 1] * work[j - 1];
      }
      tol *= tau[c_j];
      if (tol != 0.0) {
        work[c_j] = work[c_j] - tol;
        minmn = c_j + 2;
        maxmn = (((((b_A.size(0) - c_j) - 1) / 2) << 1) + c_j) + 2;
        m = maxmn - 2;
        for (int j{minmn}; j <= m; j += 2) {
          __m128d r;
          __m128d r1;
          r = _mm_loadu_pd(&b_A[(j + b_A.size(0) * c_j) - 1]);
          r1 = _mm_loadu_pd(&work[j - 1]);
          _mm_storeu_pd(&work[j - 1],
                        _mm_sub_pd(r1, _mm_mul_pd(r, _mm_set1_pd(tol))));
        }
        for (int j{maxmn}; j <= b_m; j++) {
          work[j - 1] = work[j - 1] - b_A[(j + b_A.size(0) * c_j) - 1] * tol;
        }
      }
    }
  }
  for (int j{0}; j < rankA; j++) {
    Y[jpvt[j] - 1] = work[j];
  }
  for (int j{rankA}; j >= 1; j--) {
    minmn = jpvt[j - 1];
    maxmn = b_A.size(0) * (j - 1);
    Y[minmn - 1] = Y[minmn - 1] / b_A[(j + maxmn) - 1];
    for (int c_j{0}; c_j <= j - 2; c_j++) {
      Y[jpvt[c_j] - 1] = Y[jpvt[c_j] - 1] - Y[minmn - 1] * b_A[c_j + maxmn];
    }
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for qrsolve.cpp
//
// [EOF]
//
