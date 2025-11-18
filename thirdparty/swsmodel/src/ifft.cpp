//
// File: ifft.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "ifft.h"
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "coder_array.h"
#include "omp.h"
#include <algorithm>
#include <cmath>

// Function Definitions
//
// Arguments    : const array<double, 1U> &x
//                array<creal_T, 1U> &y
// Return Type  : void
//
namespace coder {
void ifft(const array<double, 1U> &x, array<creal_T, 1U> &y)
{
  array<creal_T, 1U> b_fv;
  array<creal_T, 1U> fv;
  array<creal_T, 1U> wwc;
  array<double, 2U> costab;
  array<double, 2U> costab1q;
  array<double, 2U> sintab;
  array<double, 2U> sintabinv;
  creal_T dc;
  creal_T nt;
  double ar_tmp;
  double b_ar;
  double b_ar_tmp;
  double c_ar_tmp;
  double c_re_tmp;
  double d_ar_tmp;
  double e_re_tmp;
  double f_re_tmp;
  double g_re_tmp;
  int minNrowsNx;
  int nd2;
  minNrowsNx = x.size(0);
  if (x.size(0) == 0) {
    y.set_size(0);
  } else {
    double e;
    int N2blue;
    int b_n;
    int n;
    boolean_T useRadix2;
    useRadix2 =
        (static_cast<int>(static_cast<unsigned int>(x.size(0)) &
                          static_cast<unsigned int>(x.size(0) - 1)) == 0);
    N2blue = internal::fft::FFTImplementationCallback::get_algo_sizes(
        x.size(0), useRadix2, nd2);
    e = 6.2831853071795862 / static_cast<double>(nd2);
    n = nd2 / 2 / 2;
    costab1q.set_size(1, n + 1);
    costab1q[0] = 1.0;
    nd2 = static_cast<int>(static_cast<unsigned int>(n) >> 1) - 1;
    if (static_cast<int>(nd2 + 1 < 1600)) {
      for (int k{0}; k <= nd2; k++) {
        costab1q[k + 1] = std::cos(e * (static_cast<double>(k) + 1.0));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (int k = 0; k <= nd2; k++) {
        costab1q[k + 1] = std::cos(e * (static_cast<double>(k) + 1.0));
      }
    }
    b_n = nd2 + 2;
    if (static_cast<int>((n - nd2) - 2 < 1600)) {
      for (int b_k{b_n}; b_k < n; b_k++) {
        costab1q[b_k] = std::sin(e * static_cast<double>(n - b_k));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (int b_k = b_n; b_k < n; b_k++) {
        costab1q[b_k] = std::sin(e * static_cast<double>(n - b_k));
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
      b_n = (costab1q.size(1) - 1 < 1600);
      if (b_n) {
        for (int d_k{0}; d_k < n; d_k++) {
          sintabinv[d_k + 1] = costab1q[(n - d_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int d_k = 0; d_k < n; d_k++) {
          sintabinv[d_k + 1] = costab1q[(n - d_k) - 1];
        }
      }
      nd2 = costab1q.size(1);
      for (int f_k{nd2}; f_k <= n2; f_k++) {
        sintabinv[f_k] = costab1q[f_k - n];
      }
      if (b_n) {
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
      if ((x.size(0) != 1) &&
          ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U)) {
        int nInt2m1;
        n2 = static_cast<int>(static_cast<unsigned int>(x.size(0)) >> 1);
        nInt2m1 = (n2 + n2) - 1;
        wwc.set_size(nInt2m1);
        nd2 = 0;
        wwc[n2 - 1].re = 1.0;
        wwc[n2 - 1].im = 0.0;
        n = n2 << 1;
        for (int f_k{0}; f_k <= n2 - 2; f_k++) {
          b_n = ((f_k + 1) << 1) - 1;
          if (n - nd2 <= b_n) {
            nd2 += b_n - n;
          } else {
            nd2 += b_n;
          }
          e = 3.1415926535897931 * static_cast<double>(nd2) /
              static_cast<double>(n2);
          nt.re = std::cos(e);
          nt.im = std::sin(e);
          b_n = (n2 - f_k) - 2;
          wwc[b_n].re = nt.re;
          wwc[b_n].im = -nt.im;
        }
        b_n = nInt2m1 - 1;
        for (int f_k{b_n}; f_k >= n2; f_k--) {
          wwc[f_k] = wwc[(nInt2m1 - f_k) - 1];
        }
      } else {
        n2 = (x.size(0) + x.size(0)) - 1;
        wwc.set_size(n2);
        n = 0;
        wwc[x.size(0) - 1].re = 1.0;
        wwc[x.size(0) - 1].im = 0.0;
        nd2 = x.size(0) << 1;
        for (int f_k{0}; f_k <= minNrowsNx - 2; f_k++) {
          b_n = ((f_k + 1) << 1) - 1;
          if (nd2 - n <= b_n) {
            n += b_n - nd2;
          } else {
            n += b_n;
          }
          e = 3.1415926535897931 * static_cast<double>(n) /
              static_cast<double>(minNrowsNx);
          nt.re = std::cos(e);
          nt.im = std::sin(e);
          b_n = (x.size(0) - f_k) - 2;
          wwc[b_n].re = nt.re;
          wwc[b_n].im = -nt.im;
        }
        nd2 = n2 - 1;
        for (int f_k{nd2}; f_k >= minNrowsNx; f_k--) {
          wwc[f_k] = wwc[(n2 - f_k) - 1];
        }
      }
      y.set_size(minNrowsNx);
      if ((N2blue != 1) &&
          ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U)) {
        internal::fft::FFTImplementationCallback::b_doHalfLengthBluestein(
            x, y, x.size(0), x.size(0), N2blue, wwc, costab, sintab, costab,
            sintabinv);
      } else {
        double b_re_tmp;
        double d_re_tmp;
        double re_tmp;
        if (static_cast<int>(x.size(0) < 1600)) {
          for (int i_k{0}; i_k < minNrowsNx; i_k++) {
            nt = wwc[(minNrowsNx + i_k) - 1];
            y[i_k].re = nt.re * x[i_k];
            y[i_k].im = nt.im * -x[i_k];
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(nt)

          for (int i_k = 0; i_k < minNrowsNx; i_k++) {
            nt = wwc[(minNrowsNx + i_k) - 1];
            y[i_k].re = nt.re * x[i_k];
            y[i_k].im = nt.im * -x[i_k];
          }
        }
        nd2 = x.size(0) + 1;
        for (int f_k{nd2}; f_k <= minNrowsNx; f_k++) {
          y[f_k - 1].re = 0.0;
          y[f_k - 1].im = 0.0;
        }
        internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
            y, N2blue, costab, sintab, fv);
        internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
            wwc, N2blue, costab, sintab, b_fv);
        nd2 = fv.size(0);
        b_fv.set_size(fv.size(0));
        if (static_cast<int>(fv.size(0) < 1600)) {
          for (int i1{0}; i1 < nd2; i1++) {
            e = fv[i1].re;
            re_tmp = b_fv[i1].im;
            b_re_tmp = fv[i1].im;
            d_re_tmp = b_fv[i1].re;
            b_fv[i1].re = e * d_re_tmp - b_re_tmp * re_tmp;
            b_fv[i1].im = e * re_tmp + b_re_tmp * d_re_tmp;
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        c_re_tmp, e_re_tmp, f_re_tmp, g_re_tmp)

          for (int i1 = 0; i1 < nd2; i1++) {
            c_re_tmp = fv[i1].re;
            e_re_tmp = b_fv[i1].im;
            f_re_tmp = fv[i1].im;
            g_re_tmp = b_fv[i1].re;
            b_fv[i1].re = c_re_tmp * g_re_tmp - f_re_tmp * e_re_tmp;
            b_fv[i1].im = c_re_tmp * e_re_tmp + f_re_tmp * g_re_tmp;
          }
        }
        internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
            b_fv, N2blue, costab, sintabinv, fv);
        if (fv.size(0) > 1) {
          e = 1.0 / static_cast<double>(fv.size(0));
          nd2 = fv.size(0);
          if (static_cast<int>(fv.size(0) < 1600)) {
            for (int i2{0}; i2 < nd2; i2++) {
              fv[i2].re = e * fv[i2].re;
              fv[i2].im = e * fv[i2].im;
            }
          } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

            for (int i2 = 0; i2 < nd2; i2++) {
              fv[i2].re = e * fv[i2].re;
              fv[i2].im = e * fv[i2].im;
            }
          }
        }
        b_n = wwc.size(0);
        if (static_cast<int>((wwc.size(0) - x.size(0)) + 1 < 1600)) {
          for (int j_k{minNrowsNx}; j_k <= b_n; j_k++) {
            double ar;
            e = wwc[j_k - 1].re;
            re_tmp = fv[j_k - 1].im;
            b_re_tmp = wwc[j_k - 1].im;
            d_re_tmp = fv[j_k - 1].re;
            ar = e * d_re_tmp + b_re_tmp * re_tmp;
            e = e * re_tmp - b_re_tmp * d_re_tmp;
            if (e == 0.0) {
              nd2 = j_k - minNrowsNx;
              y[nd2].re = ar / static_cast<double>(minNrowsNx);
              y[nd2].im = 0.0;
            } else if (ar == 0.0) {
              nd2 = j_k - minNrowsNx;
              y[nd2].re = 0.0;
              y[nd2].im = e / static_cast<double>(minNrowsNx);
            } else {
              nd2 = j_k - minNrowsNx;
              y[nd2].re = ar / static_cast<double>(minNrowsNx);
              y[nd2].im = e / static_cast<double>(minNrowsNx);
            }
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        dc, ar_tmp, b_ar_tmp, c_ar_tmp, d_ar_tmp, b_ar)

          for (int j_k = minNrowsNx; j_k <= b_n; j_k++) {
            ar_tmp = wwc[j_k - 1].re;
            b_ar_tmp = fv[j_k - 1].im;
            c_ar_tmp = wwc[j_k - 1].im;
            d_ar_tmp = fv[j_k - 1].re;
            b_ar = ar_tmp * d_ar_tmp + c_ar_tmp * b_ar_tmp;
            ar_tmp = ar_tmp * b_ar_tmp - c_ar_tmp * d_ar_tmp;
            if (ar_tmp == 0.0) {
              dc.re = b_ar / static_cast<double>(minNrowsNx);
              dc.im = 0.0;
            } else if (b_ar == 0.0) {
              dc.re = 0.0;
              dc.im = ar_tmp / static_cast<double>(minNrowsNx);
            } else {
              dc.re = b_ar / static_cast<double>(minNrowsNx);
              dc.im = ar_tmp / static_cast<double>(minNrowsNx);
            }
            y[j_k - minNrowsNx] = dc;
          }
        }
      }
    } else {
      b_n = costab1q.size(1) - 1;
      n = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, n + 1);
      sintab.set_size(1, n + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      if (static_cast<int>(costab1q.size(1) - 1 < 1600)) {
        if (b_n - 1 >= 0) {
          std::copy(&costab1q[1], &costab1q[1 + b_n], &costab[1]);
        }
        for (int c_k{0}; c_k < b_n; c_k++) {
          sintab[c_k + 1] = costab1q[(b_n - c_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int c_k = 0; c_k < b_n; c_k++) {
          costab[c_k + 1] = costab1q[c_k + 1];
          sintab[c_k + 1] = costab1q[(b_n - c_k) - 1];
        }
      }
      nd2 = costab1q.size(1);
      if (static_cast<int>((n - costab1q.size(1)) + 1 < 1600)) {
        for (int e_k{nd2}; e_k <= n; e_k++) {
          costab[e_k] = -costab1q[n - e_k];
          sintab[e_k] = costab1q[e_k - b_n];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (int e_k = nd2; e_k <= n; e_k++) {
          costab[e_k] = -costab1q[n - e_k];
          sintab[e_k] = costab1q[e_k - b_n];
        }
      }
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          x, x.size(0), costab, sintab, y);
      if (y.size(0) > 1) {
        e = 1.0 / static_cast<double>(y.size(0));
        nd2 = y.size(0);
        if (static_cast<int>(y.size(0) < 1600)) {
          for (int i{0}; i < nd2; i++) {
            y[i].re = e * y[i].re;
            y[i].im = e * y[i].im;
          }
        } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

          for (int i = 0; i < nd2; i++) {
            y[i].re = e * y[i].re;
            y[i].im = e * y[i].im;
          }
        }
      }
    }
  }
}

} // namespace coder

//
// File trailer for ifft.cpp
//
// [EOF]
//
