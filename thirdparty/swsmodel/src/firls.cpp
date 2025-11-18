//
// File: firls.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "firls.h"
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "swsmodel_rtwutil.h"
#include "coder_array.h"
#include "omp.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <emmintrin.h>

// Function Declarations
static void binary_expand_op_1(coder::array<double, 1U> &in1, double in2,
                               const coder::array<double, 1U> &in3, double in4,
                               const coder::array<double, 1U> &in5);

static void binary_expand_op_2(coder::array<double, 1U> &in1, double in2,
                               const coder::array<double, 1U> &in3,
                               const coder::array<double, 1U> &in4,
                               const coder::array<double, 1U> &in5);

// Function Definitions
//
// Arguments    : coder::array<double, 1U> &in1
//                double in2
//                const coder::array<double, 1U> &in3
//                double in4
//                const coder::array<double, 1U> &in5
// Return Type  : void
//
static void binary_expand_op_1(coder::array<double, 1U> &in1, double in2,
                               const coder::array<double, 1U> &in3, double in4,
                               const coder::array<double, 1U> &in5)
{
  coder::array<double, 1U> b_in1;
  int b_in1_tmp;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    loop_ub = in1.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  b_in1.set_size(loop_ub);
  stride_0_0 = (in1.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  if (static_cast<int>(loop_ub < 1600)) {
    for (int i{0}; i < loop_ub; i++) {
      int in1_tmp;
      in1_tmp = i * stride_1_0;
      b_in1[i] =
          in1[i * stride_0_0] + (in2 * in3[in1_tmp] - in4 * in5[in1_tmp]);
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b_in1_tmp)

    for (int i = 0; i < loop_ub; i++) {
      b_in1_tmp = i * stride_1_0;
      b_in1[i] =
          in1[i * stride_0_0] + (in2 * in3[b_in1_tmp] - in4 * in5[b_in1_tmp]);
    }
  }
  in1.set_size(loop_ub);
  for (int i1{0}; i1 < loop_ub; i1++) {
    in1[i1] = b_in1[i1];
  }
}

//
// Arguments    : coder::array<double, 1U> &in1
//                double in2
//                const coder::array<double, 1U> &in3
//                const coder::array<double, 1U> &in4
//                const coder::array<double, 1U> &in5
// Return Type  : void
//
static void binary_expand_op_2(coder::array<double, 1U> &in1, double in2,
                               const coder::array<double, 1U> &in3,
                               const coder::array<double, 1U> &in4,
                               const coder::array<double, 1U> &in5)
{
  coder::array<double, 1U> b_in1;
  double d_in1_tmp;
  int c_in1_tmp;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in5.size(0) == 1) {
    loop_ub = in3.size(0);
  } else {
    loop_ub = in5.size(0);
  }
  if (loop_ub == 1) {
    loop_ub = in1.size(0);
  }
  b_in1.set_size(loop_ub);
  stride_0_0 = (in1.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_2_0 = (in5.size(0) != 1);
  if (static_cast<int>(loop_ub < 1600)) {
    for (int i{0}; i < loop_ub; i++) {
      double b_in1_tmp;
      int in1_tmp;
      in1_tmp = i * stride_1_0;
      b_in1_tmp = in5[i * stride_2_0];
      b_in1[i] = in1[i * stride_0_0] +
                 in2 * (in3[in1_tmp] - in4[in1_tmp]) / (b_in1_tmp * b_in1_tmp);
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        c_in1_tmp, d_in1_tmp)

    for (int i = 0; i < loop_ub; i++) {
      c_in1_tmp = i * stride_1_0;
      d_in1_tmp = in5[i * stride_2_0];
      b_in1[i] = in1[i * stride_0_0] + in2 * (in3[c_in1_tmp] - in4[c_in1_tmp]) /
                                           (d_in1_tmp * d_in1_tmp);
    }
  }
  in1.set_size(loop_ub);
  for (int i1{0}; i1 < loop_ub; i1++) {
    in1[i1] = b_in1[i1];
  }
}

//
// Arguments    : double N
//                const double freq[4]
//                array<double, 2U> &h
//                array<double, 1U> &a
// Return Type  : void
//
namespace coder {
void eFirls(double N, const double freq[4], array<double, 2U> &h,
            array<double, 1U> &a)
{
  static const signed char A[4]{1, 1, 0, 0};
  array<double, 2U> G;
  array<double, 2U> m;
  array<double, 2U> sinc2A;
  array<double, 2U> sinc3A;
  array<double, 2U> sinc4A;
  array<double, 1U> b_k;
  array<double, 1U> b_y;
  array<double, 1U> y;
  double F[4];
  double dv[2];
  double i2;
  double i2Map;
  double k1;
  double k2;
  double k3;
  double k4;
  double m_s;
  double tmpStorageLen;
  unsigned int b_i1;
  int b_i2;
  unsigned int i1Map;
  int idx;
  int md2;
  int nn;
  boolean_T Nodd;
  boolean_T exitg1;
  Nodd = !std::isnan(freq[0]);
  if (Nodd) {
    idx = 1;
  } else {
    idx = 0;
    md2 = 2;
    exitg1 = false;
    while ((!exitg1) && (md2 < 5)) {
      if (!std::isnan(freq[md2 - 1])) {
        idx = md2;
        exitg1 = true;
      } else {
        md2++;
      }
    }
  }
  if (idx == 0) {
    tmpStorageLen = freq[0];
  } else {
    tmpStorageLen = freq[idx - 1];
    idx++;
    for (int k{idx}; k < 5; k++) {
      m_s = freq[k - 1];
      if (tmpStorageLen < m_s) {
        tmpStorageLen = m_s;
      }
    }
  }
  if (!(tmpStorageLen > 1.0)) {
    if (Nodd) {
      idx = 1;
    } else {
      idx = 0;
      md2 = 2;
      exitg1 = false;
      while ((!exitg1) && (md2 < 5)) {
        if (!std::isnan(freq[md2 - 1])) {
          idx = md2;
          exitg1 = true;
        } else {
          md2++;
        }
      }
    }
    if (idx == 0) {
      tmpStorageLen = freq[0];
    } else {
      tmpStorageLen = freq[idx - 1];
      idx++;
      for (int k{idx}; k < 5; k++) {
        m_s = freq[k - 1];
        if (tmpStorageLen > m_s) {
          tmpStorageLen = m_s;
        }
      }
    }
    if (!(tmpStorageLen < 0.0)) {
      __m128d r;
      __m128d r1;
      double L;
      double b0;
      int b_loop_ub;
      int b_s;
      int b_vectorUB;
      int c_vectorUB;
      int i1Start;
      int loop_ub;
      int nG;
      int scalarLB;
      int ub_loop;
      int vectorUB;
      boolean_T need_matrix;
      N++;
      r = _mm_set1_pd(2.0);
      _mm_storeu_pd(&F[0], _mm_div_pd(_mm_loadu_pd(&freq[0]), r));
      _mm_storeu_pd(&F[2], _mm_div_pd(_mm_loadu_pd(&freq[2]), r));
      L = (N - 1.0) / 2.0;
      if (std::isnan(N) || std::isinf(N)) {
        tmpStorageLen = rtNaN;
      } else {
        tmpStorageLen = std::fmod(N, 2.0);
      }
      Nodd = (tmpStorageLen == 1.0);
      b0 = 0.0;
      if (!Nodd) {
        m.set_size(1, static_cast<int>(L) + 1);
        idx = static_cast<int>(L);
        md2 = ((static_cast<int>(L) + 1) / 2) << 1;
        vectorUB = md2 - 2;
        for (int k{0}; k <= vectorUB; k += 2) {
          dv[0] = k;
          dv[1] = k + 1;
          r1 = _mm_loadu_pd(&dv[0]);
          _mm_storeu_pd(&m[k], _mm_add_pd(r1, _mm_set1_pd(0.5)));
        }
        for (int k{md2}; k <= idx; k++) {
          m[k] = static_cast<double>(k) + 0.5;
        }
      } else {
        m.set_size(1, static_cast<int>(L) + 1);
        idx = static_cast<int>(L);
        for (int k{0}; k <= idx; k++) {
          m[k] = k;
        }
      }
      vectorUB = m.size(1);
      b_k.set_size(m.size(1));
      if (vectorUB - 1 >= 0) {
        std::copy(&m[0], &m[vectorUB], &b_k[0]);
      }
      need_matrix = (F[2] - F[1] != 0.0);
      nG = m.size(1);
      if (need_matrix) {
        G.set_size(m.size(1), m.size(1));
        idx = b_k.size(0) * b_k.size(0);
        if (idx - 1 >= 0) {
          std::memset(&G[0], 0,
                      static_cast<unsigned int>(idx) * sizeof(double));
        }
        tmpStorageLen = 2.0 * static_cast<double>(b_k.size(0)) - 1.0;
      } else {
        G.set_size(0, 0);
        tmpStorageLen = 0.0;
      }
      ub_loop = static_cast<int>(tmpStorageLen);
      m.set_size(1, static_cast<int>(tmpStorageLen));
      sinc2A.set_size(1, static_cast<int>(tmpStorageLen));
      sinc3A.set_size(1, static_cast<int>(tmpStorageLen));
      sinc4A.set_size(1, static_cast<int>(tmpStorageLen));
      if (Nodd) {
        i1Start = -1;
        if (b_k.size(0) < 2) {
          md2 = 0;
          vectorUB = 0;
        } else {
          md2 = 1;
        }
        idx = vectorUB - md2;
        for (int k{0}; k < idx; k++) {
          b_k[k] = b_k[md2 + k];
        }
        b_k.set_size(idx);
      } else {
        i1Start = 0;
      }
      idx = b_k.size(0);
      a.set_size(b_k.size(0));
      for (int k{0}; k < idx; k++) {
        a[k] = 0.0;
      }
      loop_ub = b_k.size(0);
      scalarLB = (b_k.size(0) / 2) << 1;
      b_vectorUB = scalarLB - 2;
      b_loop_ub = b_k.size(0);
      c_vectorUB = scalarLB - 2;
      for (int s{0}; s < 2; s++) {
        __m128d r2;
        __m128d r3;
        double b1;
        double b1_tmp;
        double m_s_tmp_tmp;
        signed char i;
        b_s = s << 1;
        m_s_tmp_tmp = F[b_s + 1];
        tmpStorageLen = F[b_s];
        N = m_s_tmp_tmp - tmpStorageLen;
        i = A[b_s];
        m_s = static_cast<double>(A[b_s + 1] - i) / N;
        b1_tmp = m_s * tmpStorageLen;
        b1 = static_cast<double>(i) - b1_tmp;
        if (Nodd) {
          b0 += b1 * N +
                m_s / 2.0 *
                    (m_s_tmp_tmp * m_s_tmp_tmp - tmpStorageLen * tmpStorageLen);
        }
        y.set_size(b_k.size(0));
        for (int k{0}; k <= b_vectorUB; k += 2) {
          r1 = _mm_loadu_pd(&b_k[k]);
          _mm_storeu_pd(&y[k], _mm_mul_pd(_mm_set1_pd(6.2831853071795862), r1));
        }
        for (int k{scalarLB}; k < loop_ub; k++) {
          y[k] = 6.2831853071795862 * b_k[k];
        }
        md2 = y.size(0);
        b_y.set_size(y.size(0));
        idx = (y.size(0) / 2) << 1;
        vectorUB = idx - 2;
        for (int k{0}; k <= vectorUB; k += 2) {
          r1 = _mm_loadu_pd(&y[k]);
          _mm_storeu_pd(&b_y[k], _mm_mul_pd(r1, _mm_set1_pd(m_s_tmp_tmp)));
        }
        for (int k{idx}; k < md2; k++) {
          b_y[k] = y[k] * m_s_tmp_tmp;
        }
        idx = b_y.size(0);
        for (int k{0}; k < idx; k++) {
          b_y[k] = std::cos(b_y[k]);
          y[k] = y[k] * F[b_s];
        }
        for (int k{0}; k < md2; k++) {
          y[k] = std::cos(y[k]);
        }
        tmpStorageLen = m_s / 39.478417604357432;
        md2 = a.size(0);
        if (b_y.size(0) == 1) {
          idx = b_k.size(0);
        } else {
          idx = b_y.size(0);
        }
        if ((b_y.size(0) == b_k.size(0)) && (a.size(0) == idx)) {
          idx = (a.size(0) / 2) << 1;
          vectorUB = idx - 2;
          for (int k{0}; k <= vectorUB; k += 2) {
            __m128d r4;
            r1 = _mm_loadu_pd(&b_y[k]);
            r2 = _mm_loadu_pd(&y[k]);
            r3 = _mm_loadu_pd(&b_k[k]);
            r4 = _mm_loadu_pd(&a[k]);
            _mm_storeu_pd(
                &a[k],
                _mm_add_pd(r4, _mm_div_pd(_mm_mul_pd(_mm_set1_pd(tmpStorageLen),
                                                     _mm_sub_pd(r1, r2)),
                                          _mm_mul_pd(r3, r3))));
          }
          for (int k{idx}; k < md2; k++) {
            a[k] = a[k] + tmpStorageLen * (b_y[k] - y[k]) / (b_k[k] * b_k[k]);
          }
        } else {
          binary_expand_op_2(a, tmpStorageLen, b_y, y, b_k);
        }
        y.set_size(b_k.size(0));
        for (int k{0}; k <= c_vectorUB; k += 2) {
          r1 = _mm_loadu_pd(&b_k[k]);
          _mm_storeu_pd(&y[k], _mm_mul_pd(r, r1));
        }
        for (int k{scalarLB}; k < b_loop_ub; k++) {
          y[k] = 2.0 * b_k[k];
        }
        md2 = y.size(0);
        b_y.set_size(y.size(0));
        idx = (y.size(0) / 2) << 1;
        vectorUB = idx - 2;
        for (int k{0}; k <= vectorUB; k += 2) {
          r1 = _mm_loadu_pd(&y[k]);
          _mm_storeu_pd(&b_y[k], _mm_mul_pd(r1, _mm_set1_pd(m_s_tmp_tmp)));
        }
        for (int k{idx}; k < md2; k++) {
          b_y[k] = y[k] * m_s_tmp_tmp;
        }
        N = m_s_tmp_tmp * (m_s * m_s_tmp_tmp + b1);
        idx = b_y.size(0);
        for (int k{0}; k < idx; k++) {
          if (std::abs(b_y[k]) < 1.0020841800044864E-292) {
            b_y[k] = 1.0;
          } else {
            tmpStorageLen = 3.1415926535897931 * b_y[k];
            tmpStorageLen = std::sin(tmpStorageLen) / tmpStorageLen;
            b_y[k] = tmpStorageLen;
          }
          y[k] = y[k] * F[b_s];
        }
        tmpStorageLen = F[b_s] * (b1_tmp + b1);
        for (int k{0}; k < md2; k++) {
          if (std::abs(y[k]) < 1.0020841800044864E-292) {
            y[k] = 1.0;
          } else {
            m_s = 3.1415926535897931 * y[k];
            m_s = std::sin(m_s) / m_s;
            y[k] = m_s;
          }
        }
        idx = a.size(0);
        if (a.size(0) == b_y.size(0)) {
          vectorUB = (a.size(0) / 2) << 1;
          md2 = vectorUB - 2;
          for (int k{0}; k <= md2; k += 2) {
            r1 = _mm_loadu_pd(&b_y[k]);
            r2 = _mm_loadu_pd(&y[k]);
            r3 = _mm_loadu_pd(&a[k]);
            _mm_storeu_pd(
                &a[k],
                _mm_add_pd(r3, _mm_sub_pd(_mm_mul_pd(_mm_set1_pd(N), r1),
                                          _mm_mul_pd(_mm_set1_pd(tmpStorageLen),
                                                     r2))));
          }
          for (int k{vectorUB}; k < idx; k++) {
            a[k] = a[k] + (N * b_y[k] - tmpStorageLen * y[k]);
          }
        } else {
          binary_expand_op_1(a, N, b_y, tmpStorageLen, y);
        }
        if (need_matrix) {
          tmpStorageLen = 2.0 * m_s_tmp_tmp;
          N = 2.0 * F[b_s];
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        k4, k3, k2, k1, i2, b_i1)

          for (int ii = 0; ii < ub_loop; ii++) {
            b_i1 = (static_cast<unsigned int>(ii) +
                    static_cast<unsigned int>(i1Start)) +
                   1U;
            i2 = (static_cast<double>(ii) + 1.0) - static_cast<double>(nG);
            k1 = tmpStorageLen * static_cast<double>(b_i1);
            k2 = N * static_cast<double>(b_i1);
            k3 = tmpStorageLen * i2;
            k4 = N * i2;
            if (std::abs(k1) < 1.0020841800044864E-292) {
              m[ii] = 1.0;
            } else {
              i2 = 3.1415926535897931 * k1;
              m[ii] = std::sin(i2) / i2;
            }
            if (std::abs(k2) < 1.0020841800044864E-292) {
              sinc2A[ii] = 1.0;
            } else {
              i2 = 3.1415926535897931 * k2;
              sinc2A[ii] = std::sin(i2) / i2;
            }
            if (std::abs(k3) < 1.0020841800044864E-292) {
              sinc3A[ii] = 1.0;
            } else {
              i2 = 3.1415926535897931 * k3;
              sinc3A[ii] = std::sin(i2) / i2;
            }
            if (std::abs(k4) < 1.0020841800044864E-292) {
              sinc4A[ii] = 1.0;
            } else {
              i2 = 3.1415926535897931 * k4;
              sinc4A[ii] = std::sin(i2) / i2;
            }
          }
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        i2Map, i1Map, nn, b_i2)

          for (int mm = 0; mm < nG; mm++) {
            for (nn = 0; nn < nG; nn++) {
              i1Map = (static_cast<unsigned int>(nn) +
                       static_cast<unsigned int>(mm)) +
                      1U;
              i2Map = static_cast<double>(nn - mm) + static_cast<double>(nG);
              b_i2 = nn + G.size(0) * mm;
              G[b_i2] =
                  G[b_i2] + 0.25 * (tmpStorageLen *
                                        (m[static_cast<int>(i1Map) - 1] +
                                         sinc3A[static_cast<int>(i2Map) - 1]) -
                                    N * (sinc2A[static_cast<int>(i1Map) - 1] +
                                         sinc4A[static_cast<int>(i2Map) - 1]));
            }
          }
        }
      }
      if (Nodd) {
        md2 = a.size(0) + 1;
        b_k.set_size(a.size(0) + 1);
        b_k[0] = b0;
        idx = a.size(0);
        for (int k{0}; k < idx; k++) {
          b_k[k + 1] = a[k];
        }
        a.set_size(md2);
        for (int k{0}; k < md2; k++) {
          a[k] = b_k[k];
        }
      }
      if (need_matrix) {
        mldivide(G, a);
      } else {
        idx = a.size(0);
        md2 = (a.size(0) / 2) << 1;
        vectorUB = md2 - 2;
        for (int k{0}; k <= vectorUB; k += 2) {
          r1 = _mm_loadu_pd(&a[k]);
          _mm_storeu_pd(&a[k], _mm_mul_pd(_mm_set1_pd(4.0), r1));
        }
        for (int k{md2}; k < idx; k++) {
          a[k] = 4.0 * a[k];
        }
        if (Nodd) {
          a[0] = a[0] / 2.0;
        }
      }
      if (Nodd) {
        if (L + 1.0 < 2.0) {
          vectorUB = 0;
          md2 = 1;
          idx = -1;
          loop_ub = 0;
          b_loop_ub = 0;
        } else {
          vectorUB = static_cast<int>(L + 1.0) - 1;
          md2 = -1;
          idx = 1;
          loop_ub = 1;
          b_loop_ub = static_cast<int>(L + 1.0);
        }
        b_s = div_s32(idx - vectorUB, md2);
        h.set_size(1, ((b_s + b_loop_ub) - loop_ub) + 2);
        if (static_cast<int>(b_s + 1 < 1600)) {
          for (int i1{0}; i1 <= b_s; i1++) {
            h[i1] = a[vectorUB + md2 * i1] / 2.0;
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

          for (int i1 = 0; i1 <= b_s; i1++) {
            h[i1] = a[vectorUB + md2 * i1] / 2.0;
          }
        }
        h[b_s + 1] = a[0];
        idx = b_loop_ub - loop_ub;
        for (int k{0}; k < idx; k++) {
          h[(k + b_s) + 2] = a[loop_ub + k] / 2.0;
        }
      } else {
        b_loop_ub = a.size(0);
        b_y.set_size(a.size(0));
        for (int k{0}; k < b_loop_ub; k++) {
          b_y[k] = a[k];
        }
        idx = a.size(0) - 1;
        md2 = a.size(0) >> 1;
        for (int k{0}; k < md2; k++) {
          tmpStorageLen = b_y[k];
          vectorUB = idx - k;
          b_y[k] = b_y[vectorUB];
          b_y[vectorUB] = tmpStorageLen;
        }
        h.set_size(1, b_y.size(0) + a.size(0));
        vectorUB = (b_y.size(0) / 2) << 1;
        idx = vectorUB - 2;
        for (int k{0}; k <= idx; k += 2) {
          r1 = _mm_loadu_pd(&b_y[k]);
          _mm_storeu_pd(&h[k], _mm_mul_pd(_mm_set1_pd(0.5), r1));
        }
        for (int k{vectorUB}; k < b_loop_ub; k++) {
          h[k] = 0.5 * b_y[k];
        }
        idx = vectorUB - 2;
        for (int k{0}; k <= idx; k += 2) {
          r1 = _mm_loadu_pd(&a[k]);
          _mm_storeu_pd(&h[k + b_y.size(0)], _mm_mul_pd(_mm_set1_pd(0.5), r1));
        }
        for (int k{vectorUB}; k < b_loop_ub; k++) {
          h[k + b_y.size(0)] = 0.5 * a[k];
        }
      }
      a.set_size(1);
      a[0] = 1.0;
    }
  }
}

} // namespace coder

//
// File trailer for firls.cpp
//
// [EOF]
//
