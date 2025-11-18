//
// File: fft.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "fft.h"
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const array<double, 1U> &x
//                double varargin_1
//                array<creal_T, 1U> &y
// Return Type  : void
//
namespace coder {
void fft(const array<double, 1U> &x, double varargin_1, array<creal_T, 1U> &y)
{
  array<creal_T, 1U> b_fv;
  array<creal_T, 1U> fv;
  array<creal_T, 1U> wwc;
  array<double, 2U> costab;
  array<double, 2U> costab1q;
  array<double, 2U> sintab;
  array<double, 2U> sintabinv;
  creal_T nt;
  double c_re_tmp;
  double d;
  double d1;
  double d2;
  double d3;
  double e_re_tmp;
  double f_re_tmp;
  double g_re_tmp;
  int i2;
  int nd2;
  int nfft;
  nfft = static_cast<int>(varargin_1);
  if ((x.size(0) == 0) || (static_cast<int>(varargin_1) == 0)) {
    y.set_size(static_cast<int>(varargin_1));
    for (int k{0}; k < nfft; k++) {
      y[k].re = 0.0;
      y[k].im = 0.0;
    }
  } else {
    double e;
    int N2blue;
    int minNrowsNx;
    int n;
    boolean_T useRadix2;
    useRadix2 = ((static_cast<int>(varargin_1) > 0) &&
                 ((static_cast<int>(varargin_1) &
                   (static_cast<int>(varargin_1) - 1)) == 0));
    N2blue = internal::fft::FFTImplementationCallback::get_algo_sizes(
        static_cast<int>(varargin_1), useRadix2, nd2);
    e = 6.2831853071795862 / static_cast<double>(nd2);
    n = nd2 / 2 / 2;
    costab1q.set_size(1, n + 1);
    costab1q[0] = 1.0;
    nd2 = static_cast<int>(static_cast<unsigned int>(n) >> 1) - 1;
    if (static_cast<int>(nd2 + 1 < 1600)) {
      for (int b_k{0}; b_k <= nd2; b_k++) {
        costab1q[b_k + 1] = std::cos(e * (static_cast<double>(b_k) + 1.0));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (int b_k = 0; b_k <= nd2; b_k++) {
        costab1q[b_k + 1] = std::cos(e * (static_cast<double>(b_k) + 1.0));
      }
    }
    minNrowsNx = nd2 + 2;
    if (static_cast<int>((n - nd2) - 2 < 1600)) {
      for (int c_k{minNrowsNx}; c_k < n; c_k++) {
        costab1q[c_k] = std::sin(e * static_cast<double>(n - c_k));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (int c_k = minNrowsNx; c_k < n; c_k++) {
        costab1q[c_k] = std::sin(e * static_cast<double>(n - c_k));
      }
    }
    costab1q[n] = 0.0;
    if (!useRadix2) {
      int n2;
      n = costab1q.size(1) - 1;
      n2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, n2 + 1);
      sintab.set_size(1, n2 + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, n2 + 1);
      minNrowsNx = (costab1q.size(1) - 1 < 1600);
      if (minNrowsNx) {
        for (int e_k{0}; e_k < n; e_k++) {
          sintabinv[e_k + 1] = costab1q[(n - e_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int e_k = 0; e_k < n; e_k++) {
          sintabinv[e_k + 1] = costab1q[(n - e_k) - 1];
        }
      }
      nd2 = costab1q.size(1);
      for (int k{nd2}; k <= n2; k++) {
        sintabinv[k] = costab1q[k - n];
      }
      if (minNrowsNx) {
        for (int g_k{0}; g_k < n; g_k++) {
          costab[g_k + 1] = costab1q[g_k + 1];
          sintab[g_k + 1] = -costab1q[(n - g_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int g_k = 0; g_k < n; g_k++) {
          costab[g_k + 1] = costab1q[g_k + 1];
          sintab[g_k + 1] = -costab1q[(n - g_k) - 1];
        }
      }
      if (static_cast<int>((n2 - costab1q.size(1)) + 1 < 1600)) {
        for (int h_k{nd2}; h_k <= n2; h_k++) {
          costab[h_k] = -costab1q[n2 - h_k];
          sintab[h_k] = -costab1q[h_k - n];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int h_k = nd2; h_k <= n2; h_k++) {
          costab[h_k] = -costab1q[n2 - h_k];
          sintab[h_k] = -costab1q[h_k - n];
        }
      }
      if ((static_cast<int>(varargin_1) != 1) &&
          ((static_cast<unsigned int>(static_cast<int>(varargin_1)) & 1U) ==
           0U)) {
        int nInt2m1;
        n2 = static_cast<int>(varargin_1) / 2;
        nInt2m1 = (n2 + n2) - 1;
        wwc.set_size(nInt2m1);
        nd2 = 0;
        wwc[n2 - 1].re = 1.0;
        wwc[n2 - 1].im = 0.0;
        n = n2 << 1;
        for (int k{0}; k <= n2 - 2; k++) {
          minNrowsNx = ((k + 1) << 1) - 1;
          if (n - nd2 <= minNrowsNx) {
            nd2 += minNrowsNx - n;
          } else {
            nd2 += minNrowsNx;
          }
          e = -3.1415926535897931 * static_cast<double>(nd2) /
              static_cast<double>(n2);
          nt.re = std::cos(e);
          nt.im = std::sin(e);
          minNrowsNx = (n2 - k) - 2;
          wwc[minNrowsNx].re = nt.re;
          wwc[minNrowsNx].im = -nt.im;
        }
        minNrowsNx = nInt2m1 - 1;
        for (int k{minNrowsNx}; k >= n2; k--) {
          wwc[k] = wwc[(nInt2m1 - k) - 1];
        }
      } else {
        n2 = (static_cast<int>(varargin_1) + static_cast<int>(varargin_1)) - 1;
        wwc.set_size(n2);
        n = 0;
        wwc[static_cast<int>(varargin_1) - 1].re = 1.0;
        wwc[static_cast<int>(varargin_1) - 1].im = 0.0;
        nd2 = static_cast<int>(varargin_1) << 1;
        for (int k{0}; k <= nfft - 2; k++) {
          minNrowsNx = ((k + 1) << 1) - 1;
          if (nd2 - n <= minNrowsNx) {
            n += minNrowsNx - nd2;
          } else {
            n += minNrowsNx;
          }
          e = -3.1415926535897931 * static_cast<double>(n) /
              static_cast<double>(static_cast<int>(varargin_1));
          nt.re = std::cos(e);
          nt.im = std::sin(e);
          minNrowsNx = (static_cast<int>(varargin_1) - k) - 2;
          wwc[minNrowsNx].re = nt.re;
          wwc[minNrowsNx].im = -nt.im;
        }
        nd2 = n2 - 1;
        for (int k{nd2}; k >= nfft; k--) {
          wwc[k] = wwc[(n2 - k) - 1];
        }
      }
      y.set_size(nfft);
      if (static_cast<int>(varargin_1) > x.size(0)) {
        y.set_size(nfft);
        for (int k{0}; k < nfft; k++) {
          y[k].re = 0.0;
          y[k].im = 0.0;
        }
      }
      if ((N2blue != 1) &&
          ((static_cast<unsigned int>(static_cast<int>(varargin_1)) & 1U) ==
           0U)) {
        internal::fft::FFTImplementationCallback::doHalfLengthBluestein(
            x, y, x.size(0), static_cast<int>(varargin_1), N2blue, wwc, costab,
            sintab, costab, sintabinv);
      } else {
        double b_re_tmp;
        double d_re_tmp;
        double re_tmp;
        nd2 = static_cast<int>(varargin_1);
        minNrowsNx = x.size(0);
        if (nd2 <= minNrowsNx) {
          minNrowsNx = nd2;
        }
        if (static_cast<int>(minNrowsNx < 1600)) {
          for (int i_k{0}; i_k < minNrowsNx; i_k++) {
            nt = wwc[(static_cast<int>(varargin_1) + i_k) - 1];
            y[i_k].re = nt.re * x[i_k];
            y[i_k].im = nt.im * -x[i_k];
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(nt)

          for (int i_k = 0; i_k < minNrowsNx; i_k++) {
            nt = wwc[(static_cast<int>(varargin_1) + i_k) - 1];
            y[i_k].re = nt.re * x[i_k];
            y[i_k].im = nt.im * -x[i_k];
          }
        }
        nd2 = minNrowsNx + 1;
        for (int k{nd2}; k <= nfft; k++) {
          y[k - 1].re = 0.0;
          y[k - 1].im = 0.0;
        }
        internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
            y, N2blue, costab, sintab, fv);
        internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
            wwc, N2blue, costab, sintab, b_fv);
        nd2 = fv.size(0);
        b_fv.set_size(fv.size(0));
        if (static_cast<int>(fv.size(0) < 1600)) {
          for (int i{0}; i < nd2; i++) {
            e = fv[i].re;
            re_tmp = b_fv[i].im;
            b_re_tmp = fv[i].im;
            d_re_tmp = b_fv[i].re;
            b_fv[i].re = e * d_re_tmp - b_re_tmp * re_tmp;
            b_fv[i].im = e * re_tmp + b_re_tmp * d_re_tmp;
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        c_re_tmp, e_re_tmp, f_re_tmp, g_re_tmp)

          for (int i = 0; i < nd2; i++) {
            c_re_tmp = fv[i].re;
            e_re_tmp = b_fv[i].im;
            f_re_tmp = fv[i].im;
            g_re_tmp = b_fv[i].re;
            b_fv[i].re = c_re_tmp * g_re_tmp - f_re_tmp * e_re_tmp;
            b_fv[i].im = c_re_tmp * e_re_tmp + f_re_tmp * g_re_tmp;
          }
        }
        internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
            b_fv, N2blue, costab, sintabinv, fv);
        if (fv.size(0) > 1) {
          e = 1.0 / static_cast<double>(fv.size(0));
          nd2 = fv.size(0);
          if (static_cast<int>(fv.size(0) < 1600)) {
            for (int i1{0}; i1 < nd2; i1++) {
              fv[i1].re = e * fv[i1].re;
              fv[i1].im = e * fv[i1].im;
            }
          } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

            for (int i1 = 0; i1 < nd2; i1++) {
              fv[i1].re = e * fv[i1].re;
              fv[i1].im = e * fv[i1].im;
            }
          }
        }
        minNrowsNx = wwc.size(0);
        if (static_cast<int>((wwc.size(0) - static_cast<int>(varargin_1)) + 1 <
                             1600)) {
          for (int j_k{nfft}; j_k <= minNrowsNx; j_k++) {
            e = wwc[j_k - 1].re;
            re_tmp = fv[j_k - 1].im;
            b_re_tmp = wwc[j_k - 1].im;
            d_re_tmp = fv[j_k - 1].re;
            nd2 = j_k - static_cast<int>(varargin_1);
            y[nd2].re = e * d_re_tmp + b_re_tmp * re_tmp;
            y[nd2].im = e * re_tmp - b_re_tmp * d_re_tmp;
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d, d1, d2, d3, i2)

          for (int j_k = nfft; j_k <= minNrowsNx; j_k++) {
            d = wwc[j_k - 1].re;
            d1 = fv[j_k - 1].im;
            d2 = wwc[j_k - 1].im;
            d3 = fv[j_k - 1].re;
            i2 = j_k - nfft;
            y[i2].re = d * d3 + d2 * d1;
            y[i2].im = d * d1 - d2 * d3;
          }
        }
      }
    } else {
      minNrowsNx = costab1q.size(1) - 1;
      n = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, n + 1);
      sintab.set_size(1, n + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      if (static_cast<int>(costab1q.size(1) - 1 < 1600)) {
        for (int d_k{0}; d_k < minNrowsNx; d_k++) {
          costab[d_k + 1] = costab1q[d_k + 1];
          sintab[d_k + 1] = -costab1q[(minNrowsNx - d_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int d_k = 0; d_k < minNrowsNx; d_k++) {
          costab[d_k + 1] = costab1q[d_k + 1];
          sintab[d_k + 1] = -costab1q[(minNrowsNx - d_k) - 1];
        }
      }
      nd2 = costab1q.size(1);
      if (static_cast<int>((n - costab1q.size(1)) + 1 < 1600)) {
        for (int f_k{nd2}; f_k <= n; f_k++) {
          costab[f_k] = -costab1q[n - f_k];
          sintab[f_k] = -costab1q[f_k - minNrowsNx];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int f_k = nd2; f_k <= n; f_k++) {
          costab[f_k] = -costab1q[n - f_k];
          sintab[f_k] = -costab1q[f_k - minNrowsNx];
        }
      }
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          x, static_cast<int>(varargin_1), costab, sintab, y);
    }
  }
}

} // namespace coder

//
// File trailer for fft.cpp
//
// [EOF]
//
