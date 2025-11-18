//
// File: FFTImplementationCallback.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const array<double, 1U> &x
//                array<creal_T, 1U> &y
//                int unsigned_nRows
//                const array<double, 2U> &costab
//                const array<double, 2U> &sintab
// Return Type  : void
//
namespace coder {
namespace internal {
namespace fft {
void FFTImplementationCallback::doHalfLengthRadix2(
    const array<double, 1U> &x, array<creal_T, 1U> &y, int unsigned_nRows,
    const array<double, 2U> &costab, const array<double, 2U> &sintab)
{
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<double, 2U> hcostab;
  array<double, 2U> hsintab;
  array<int, 2U> wrapIndex;
  array<int, 1U> bitrevIndex;
  double b_im;
  double b_re;
  double im;
  double re;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_re;
  int b_i;
  int hszCostab;
  int iDelta;
  int iDelta2;
  int iheight;
  int iy;
  int ju;
  int k;
  int nRows;
  int nRowsD2;
  int u0;
  boolean_T tst;
  nRows = unsigned_nRows / 2;
  u0 = y.size(0);
  if (u0 > nRows) {
    u0 = nRows;
  }
  iDelta = u0 - 2;
  iDelta2 = nRows - 2;
  nRowsD2 = nRows / 2;
  k = nRowsD2 / 2;
  hszCostab = static_cast<int>(static_cast<unsigned int>(costab.size(1)) >> 1);
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  if (static_cast<int>(hszCostab < 1600)) {
    for (int i{0}; i < hszCostab; i++) {
      b_i = ((i + 1) << 1) - 2;
      hcostab[i] = costab[b_i];
      hsintab[i] = sintab[b_i];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b_i)

    for (int i = 0; i < hszCostab; i++) {
      b_i = ((i + 1) << 1) - 2;
      hcostab[i] = costab[b_i];
      hsintab[i] = sintab[b_i];
    }
  }
  reconVar1.set_size(nRows);
  reconVar2.set_size(nRows);
  wrapIndex.set_size(1, nRows);
  ju = 0;
  iy = 1;
  bitrevIndex.set_size(nRows);
  for (int c_i{0}; c_i < nRows; c_i++) {
    re = sintab[c_i];
    im = costab[c_i];
    reconVar1[c_i].re = re + 1.0;
    reconVar1[c_i].im = -im;
    reconVar2[c_i].re = 1.0 - re;
    reconVar2[c_i].im = im;
    if (c_i != 0) {
      wrapIndex[c_i] = (nRows - c_i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
    bitrevIndex[c_i] = 0;
  }
  for (int c_i{0}; c_i <= iDelta; c_i++) {
    bitrevIndex[c_i] = iy;
    iy = nRows;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }
    iy = ju + 1;
  }
  bitrevIndex[u0 - 1] = iy;
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    tst = true;
    hszCostab = x.size(0);
  } else if (x.size(0) >= unsigned_nRows) {
    tst = true;
    hszCostab = unsigned_nRows;
  } else {
    tst = false;
    hszCostab = x.size(0) - 1;
  }
  if (hszCostab > unsigned_nRows) {
    hszCostab = unsigned_nRows;
  }
  ju = hszCostab / 2;
  for (int c_i{0}; c_i < ju; c_i++) {
    hszCostab = c_i << 1;
    y[bitrevIndex[c_i] - 1].re = x[hszCostab];
    y[bitrevIndex[c_i] - 1].im = x[hszCostab + 1];
  }
  if (!tst) {
    if (ju - 1 < 0) {
      hszCostab = 0;
    } else {
      hszCostab = ju << 1;
    }
    y[bitrevIndex[ju] - 1].re = x[hszCostab];
    y[bitrevIndex[ju] - 1].im = 0.0;
  }
  if (nRows > 1) {
    for (int c_i{0}; c_i <= iDelta2; c_i += 2) {
      re = y[c_i + 1].re;
      im = y[c_i + 1].im;
      temp_re = re;
      temp_im = im;
      b_re = y[c_i].re;
      b_im = y[c_i].im;
      re = b_re - re;
      im = b_im - im;
      y[c_i + 1].re = re;
      y[c_i + 1].im = im;
      b_re += temp_re;
      b_im += temp_im;
      y[c_i].re = b_re;
      y[c_i].im = b_im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int d_i;
    for (d_i = 0; d_i < iheight; d_i += iDelta2) {
      hszCostab = d_i + iDelta;
      temp_re = y[hszCostab].re;
      temp_im = y[hszCostab].im;
      y[hszCostab].re = y[d_i].re - temp_re;
      y[hszCostab].im = y[d_i].im - temp_im;
      y[d_i].re = y[d_i].re + temp_re;
      y[d_i].im = y[d_i].im + temp_im;
    }
    hszCostab = 1;
    for (iy = k; iy < nRowsD2; iy += k) {
      temp2_re = hcostab[iy];
      temp2_im = hsintab[iy];
      d_i = hszCostab;
      ju = hszCostab + iheight;
      while (d_i < ju) {
        u0 = d_i + iDelta;
        re = y[u0].im;
        im = y[u0].re;
        temp_re = temp2_re * im - temp2_im * re;
        temp_im = temp2_re * re + temp2_im * im;
        y[u0].re = y[d_i].re - temp_re;
        y[u0].im = y[d_i].im - temp_im;
        y[d_i].re = y[d_i].re + temp_re;
        y[d_i].im = y[d_i].im + temp_im;
        d_i += iDelta2;
      }
      hszCostab++;
    }
    k = static_cast<int>(static_cast<unsigned int>(k) >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
  iy = nRows / 2;
  temp_re = y[0].re;
  temp_im = y[0].im;
  im = y[0].re * reconVar1[0].re;
  b_re = y[0].re * reconVar1[0].im;
  b_im = -y[0].im;
  temp2_re = temp_re * reconVar2[0].re;
  re = temp_re * reconVar2[0].im;
  y[0].re = 0.5 * ((im - y[0].im * reconVar1[0].im) +
                   (temp2_re - b_im * reconVar2[0].im));
  y[0].im = 0.5 * ((b_re + y[0].im * reconVar1[0].re) +
                   (re + b_im * reconVar2[0].re));
  y[nRows].re = 0.5 * ((temp2_re - temp_im * reconVar2[0].im) +
                       (im - b_im * reconVar1[0].im));
  y[nRows].im = 0.5 * ((re + temp_im * reconVar2[0].re) +
                       (b_re + b_im * reconVar1[0].re));
  for (int c_i{2}; c_i <= iy; c_i++) {
    temp_re = y[c_i - 1].re;
    temp_im = y[c_i - 1].im;
    hszCostab = wrapIndex[c_i - 1];
    temp2_re = y[hszCostab - 1].re;
    temp2_im = y[hszCostab - 1].im;
    re = reconVar1[c_i - 1].im;
    im = reconVar1[c_i - 1].re;
    b_re = reconVar2[c_i - 1].im;
    b_im = reconVar2[c_i - 1].re;
    y[c_i - 1].re = 0.5 * ((temp_re * im - temp_im * re) +
                           (temp2_re * b_im - -temp2_im * b_re));
    y[c_i - 1].im = 0.5 * ((temp_re * re + temp_im * im) +
                           (temp2_re * b_re + -temp2_im * b_im));
    ju = (nRows + c_i) - 1;
    y[ju].re = 0.5 * ((temp_re * b_im - temp_im * b_re) +
                      (temp2_re * im - -temp2_im * re));
    y[ju].im = 0.5 * ((temp_re * b_re + temp_im * b_im) +
                      (temp2_re * re + -temp2_im * im));
    re = reconVar1[hszCostab - 1].im;
    im = reconVar1[hszCostab - 1].re;
    b_re = reconVar2[hszCostab - 1].im;
    b_im = reconVar2[hszCostab - 1].re;
    y[hszCostab - 1].re = 0.5 * ((temp2_re * im - temp2_im * re) +
                                 (temp_re * b_im - -temp_im * b_re));
    y[hszCostab - 1].im = 0.5 * ((temp2_re * re + temp2_im * im) +
                                 (temp_re * b_re + -temp_im * b_im));
    hszCostab = (hszCostab + nRows) - 1;
    y[hszCostab].re = 0.5 * ((temp2_re * b_im - temp2_im * b_re) +
                             (temp_re * im - -temp_im * re));
    y[hszCostab].im = 0.5 * ((temp2_re * b_re + temp2_im * b_im) +
                             (temp_re * re + -temp_im * im));
  }
  if (iy != 0) {
    double b_y_re_tmp;
    double y_re_tmp;
    temp_re = y[iy].re;
    temp_im = y[iy].im;
    im = reconVar1[iy].im;
    b_re = reconVar1[iy].re;
    b_im = temp_re * b_re;
    temp2_re = temp_re * im;
    temp2_im = reconVar2[iy].im;
    y_re_tmp = reconVar2[iy].re;
    b_y_re_tmp = temp_re * y_re_tmp;
    re = temp_re * temp2_im;
    y[iy].re =
        0.5 * ((b_im - temp_im * im) + (b_y_re_tmp - -temp_im * temp2_im));
    y[iy].im = 0.5 * ((temp2_re + temp_im * b_re) + (re + -temp_im * y_re_tmp));
    hszCostab = nRows + iy;
    y[hszCostab].re =
        0.5 * ((b_y_re_tmp - temp_im * temp2_im) + (b_im - -temp_im * im));
    y[hszCostab].im =
        0.5 * ((re + temp_im * y_re_tmp) + (temp2_re + -temp_im * b_re));
  }
}

//
// Arguments    : int nRows
//                array<double, 2U> &costab
//                array<double, 2U> &sintab
//                array<double, 2U> &sintabinv
// Return Type  : void
//
void FFTImplementationCallback::generate_twiddle_tables(
    int nRows, array<double, 2U> &costab, array<double, 2U> &sintab,
    array<double, 2U> &sintabinv)
{
  array<double, 2U> costab1q;
  double e;
  int i;
  int n;
  int n2;
  int nd2;
  e = 6.2831853071795862 / static_cast<double>(nRows);
  n = nRows / 2 / 2;
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
  i = nd2 + 2;
  if (static_cast<int>((n - nd2) - 2 < 1600)) {
    for (int b_k{i}; b_k < n; b_k++) {
      costab1q[b_k] = std::sin(e * static_cast<double>(n - b_k));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (int b_k = i; b_k < n; b_k++) {
      costab1q[b_k] = std::sin(e * static_cast<double>(n - b_k));
    }
  }
  costab1q[n] = 0.0;
  n = costab1q.size(1) - 1;
  n2 = (costab1q.size(1) - 1) << 1;
  costab.set_size(1, n2 + 1);
  sintab.set_size(1, n2 + 1);
  costab[0] = 1.0;
  sintab[0] = 0.0;
  sintabinv.set_size(1, n2 + 1);
  nd2 = (costab1q.size(1) - 1 < 1600);
  if (nd2) {
    for (int c_k{0}; c_k < n; c_k++) {
      sintabinv[c_k + 1] = costab1q[(n - c_k) - 1];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (int c_k = 0; c_k < n; c_k++) {
      sintabinv[c_k + 1] = costab1q[(n - c_k) - 1];
    }
  }
  i = costab1q.size(1);
  for (int d_k{i}; d_k <= n2; d_k++) {
    sintabinv[d_k] = costab1q[d_k - n];
  }
  if (nd2) {
    for (int e_k{0}; e_k < n; e_k++) {
      costab[e_k + 1] = costab1q[e_k + 1];
      sintab[e_k + 1] = -costab1q[(n - e_k) - 1];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (int e_k = 0; e_k < n; e_k++) {
      costab[e_k + 1] = costab1q[e_k + 1];
      sintab[e_k + 1] = -costab1q[(n - e_k) - 1];
    }
  }
  if (static_cast<int>((n2 - costab1q.size(1)) + 1 < 1600)) {
    for (int f_k{i}; f_k <= n2; f_k++) {
      costab[f_k] = -costab1q[n2 - f_k];
      sintab[f_k] = -costab1q[f_k - n];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (int f_k = i; f_k <= n2; f_k++) {
      costab[f_k] = -costab1q[n2 - f_k];
      sintab[f_k] = -costab1q[f_k - n];
    }
  }
}

//
// Arguments    : const array<double, 2U> &costab
//                const array<double, 2U> &sintab
//                const array<double, 2U> &costabinv
//                const array<double, 2U> &sintabinv
//                array<double, 2U> &hcostab
//                array<double, 2U> &hsintab
//                array<double, 2U> &hcostabinv
//                array<double, 2U> &hsintabinv
// Return Type  : void
//
void FFTImplementationCallback::get_half_twiddle_tables(
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    const array<double, 2U> &costabinv, const array<double, 2U> &sintabinv,
    array<double, 2U> &hcostab, array<double, 2U> &hsintab,
    array<double, 2U> &hcostabinv, array<double, 2U> &hsintabinv)
{
  int b_i;
  int hszCostab;
  hszCostab = static_cast<int>(static_cast<unsigned int>(costab.size(1)) >> 1);
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  hcostabinv.set_size(1, hszCostab);
  hsintabinv.set_size(1, hszCostab);
  if (static_cast<int>(hszCostab < 1600)) {
    for (int i{0}; i < hszCostab; i++) {
      b_i = ((i + 1) << 1) - 2;
      hcostab[i] = costab[b_i];
      hsintab[i] = sintab[b_i];
      hcostabinv[i] = costabinv[b_i];
      hsintabinv[i] = sintabinv[b_i];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b_i)

    for (int i = 0; i < hszCostab; i++) {
      b_i = ((i + 1) << 1) - 2;
      hcostab[i] = costab[b_i];
      hsintab[i] = sintab[b_i];
      hcostabinv[i] = costabinv[b_i];
      hsintabinv[i] = sintabinv[b_i];
    }
  }
}

//
// Arguments    : const array<double, 1U> &x
//                array<creal_T, 1U> &y
//                int nrowsx
//                int nRows
//                int nfft
//                const array<creal_T, 1U> &wwc
//                const array<double, 2U> &costab
//                const array<double, 2U> &sintab
//                const array<double, 2U> &costabinv
//                const array<double, 2U> &sintabinv
// Return Type  : void
//
void FFTImplementationCallback::b_doHalfLengthBluestein(
    const array<double, 1U> &x, array<creal_T, 1U> &y, int nrowsx, int nRows,
    int nfft, const array<creal_T, 1U> &wwc, const array<double, 2U> &costab,
    const array<double, 2U> &sintab, const array<double, 2U> &costabinv,
    const array<double, 2U> &sintabinv)
{
  array<creal_T, 1U> b_fv;
  array<creal_T, 1U> fv;
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<creal_T, 1U> ytmp;
  array<double, 2U> a__1;
  array<double, 2U> costable;
  array<double, 2U> hcostabinv;
  array<double, 2U> hsintab;
  array<double, 2U> hsintabinv;
  array<double, 2U> sintable;
  array<int, 2U> wrapIndex;
  creal_T a;
  creal_T b;
  creal_T dc;
  double ai;
  double ar;
  double ar_tmp;
  double b_ar;
  double b_ar_tmp;
  double b_re_tmp;
  double b_ytmp_im;
  double c_ar_tmp;
  double c_re_tmp;
  double c_ytmp_re_tmp;
  double d_ar_tmp;
  double d_re_tmp;
  double d_ytmp_re_tmp;
  double e_re_tmp;
  double e_ytmp_re_tmp;
  double f_re_tmp;
  double f_ytmp_re_tmp;
  double g_re_tmp;
  double g_ytmp_re_tmp;
  double h_re_tmp;
  double h_ytmp_re_tmp;
  double im;
  double re;
  double re_tmp;
  double ytmp_re;
  int b_tmp;
  int hnRows;
  int i2;
  int minHnrowsNxBy2;
  int u0;
  boolean_T nxeven;
  hnRows = nRows / 2;
  ytmp.set_size(hnRows);
  if (hnRows > nrowsx) {
    ytmp.set_size(hnRows);
    if (hnRows - 1 >= 0) {
      std::memset(&ytmp[0], 0,
                  static_cast<unsigned int>(hnRows) * sizeof(creal_T));
    }
  }
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    nxeven = true;
    u0 = x.size(0);
  } else if (x.size(0) >= nRows) {
    nxeven = true;
    u0 = nRows;
  } else {
    nxeven = false;
    u0 = x.size(0) - 1;
  }
  FFTImplementationCallback::generate_twiddle_tables(nRows << 1, costable,
                                                     sintable, a__1);
  FFTImplementationCallback::get_half_twiddle_tables(costab, sintab, costabinv,
                                                     sintabinv, a__1, hsintab,
                                                     hcostabinv, hsintabinv);
  reconVar1.set_size(hnRows);
  reconVar2.set_size(hnRows);
  wrapIndex.set_size(1, hnRows);
  for (int i{0}; i < hnRows; i++) {
    minHnrowsNxBy2 = i << 1;
    re_tmp = sintable[minHnrowsNxBy2];
    b_re_tmp = costable[minHnrowsNxBy2];
    reconVar1[i].re = 1.0 - re_tmp;
    reconVar1[i].im = -b_re_tmp;
    reconVar2[i].re = re_tmp + 1.0;
    reconVar2[i].im = b_re_tmp;
    if (i != 0) {
      wrapIndex[i] = (hnRows - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
  }
  if (u0 > nRows) {
    u0 = nRows;
  }
  minHnrowsNxBy2 = u0 / 2 - 1;
  if (static_cast<int>(minHnrowsNxBy2 + 1 < 1600)) {
    for (int k1{0}; k1 <= minHnrowsNxBy2; k1++) {
      a = wwc[(hnRows + k1) - 1];
      u0 = k1 << 1;
      b.re = x[u0];
      b.im = x[u0 + 1];
      ytmp[k1].re = a.re * b.re + a.im * b.im;
      ytmp[k1].im = a.re * b.im - a.im * b.re;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b, a, b_tmp)

    for (int k1 = 0; k1 <= minHnrowsNxBy2; k1++) {
      a = wwc[(hnRows + k1) - 1];
      b_tmp = k1 << 1;
      b.re = x[b_tmp];
      b.im = x[b_tmp + 1];
      ytmp[k1].re = a.re * b.re + a.im * b.im;
      ytmp[k1].im = a.re * b.im - a.im * b.re;
    }
  }
  if (!nxeven) {
    a = wwc[hnRows + minHnrowsNxBy2];
    if (minHnrowsNxBy2 < 0) {
      u0 = 0;
    } else {
      u0 = (minHnrowsNxBy2 + 1) << 1;
    }
    b.re = x[u0];
    ytmp[minHnrowsNxBy2 + 1].re = a.re * b.re + a.im * 0.0;
    ytmp[minHnrowsNxBy2 + 1].im = a.re * 0.0 - a.im * b.re;
    if (minHnrowsNxBy2 + 3 <= hnRows) {
      u0 = minHnrowsNxBy2 + 3;
      if (u0 <= hnRows) {
        std::memset(&ytmp[u0 + -1], 0,
                    static_cast<unsigned int>((hnRows - u0) + 1) *
                        sizeof(creal_T));
      }
    }
  } else if (minHnrowsNxBy2 + 2 <= hnRows) {
    u0 = minHnrowsNxBy2 + 2;
    if (u0 <= hnRows) {
      std::memset(&ytmp[u0 + -1], 0,
                  static_cast<unsigned int>((hnRows - u0) + 1) *
                      sizeof(creal_T));
    }
  }
  minHnrowsNxBy2 = nfft / 2;
  FFTImplementationCallback::r2br_r2dit_trig_impl(ytmp, minHnrowsNxBy2, a__1,
                                                  hsintab, fv);
  FFTImplementationCallback::r2br_r2dit_trig_impl(wwc, minHnrowsNxBy2, a__1,
                                                  hsintab, b_fv);
  u0 = fv.size(0);
  b_fv.set_size(fv.size(0));
  if (static_cast<int>(fv.size(0) < 1600)) {
    for (int b_i{0}; b_i < u0; b_i++) {
      re_tmp = fv[b_i].re;
      b_re_tmp = b_fv[b_i].im;
      c_re_tmp = fv[b_i].im;
      e_re_tmp = b_fv[b_i].re;
      b_fv[b_i].re = re_tmp * e_re_tmp - c_re_tmp * b_re_tmp;
      b_fv[b_i].im = re_tmp * b_re_tmp + c_re_tmp * e_re_tmp;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d_re_tmp, f_re_tmp, g_re_tmp, h_re_tmp)

    for (int b_i = 0; b_i < u0; b_i++) {
      d_re_tmp = fv[b_i].re;
      f_re_tmp = b_fv[b_i].im;
      g_re_tmp = fv[b_i].im;
      h_re_tmp = b_fv[b_i].re;
      b_fv[b_i].re = d_re_tmp * h_re_tmp - g_re_tmp * f_re_tmp;
      b_fv[b_i].im = d_re_tmp * f_re_tmp + g_re_tmp * h_re_tmp;
    }
  }
  FFTImplementationCallback::r2br_r2dit_trig_impl(b_fv, minHnrowsNxBy2,
                                                  hcostabinv, hsintabinv, fv);
  if (fv.size(0) > 1) {
    re_tmp = 1.0 / static_cast<double>(fv.size(0));
    u0 = fv.size(0);
    if (static_cast<int>(fv.size(0) < 1600)) {
      for (int i1{0}; i1 < u0; i1++) {
        fv[i1].re = re_tmp * fv[i1].re;
        fv[i1].im = re_tmp * fv[i1].im;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (int i1 = 0; i1 < u0; i1++) {
        fv[i1].re = re_tmp * fv[i1].re;
        fv[i1].im = re_tmp * fv[i1].im;
      }
    }
  }
  minHnrowsNxBy2 = wwc.size(0);
  if (static_cast<int>((wwc.size(0) - hnRows) + 1 < 1600)) {
    for (int k{hnRows}; k <= minHnrowsNxBy2; k++) {
      re_tmp = wwc[k - 1].re;
      b_re_tmp = fv[k - 1].im;
      c_re_tmp = wwc[k - 1].im;
      e_re_tmp = fv[k - 1].re;
      ar = re_tmp * e_re_tmp + c_re_tmp * b_re_tmp;
      re_tmp = re_tmp * b_re_tmp - c_re_tmp * e_re_tmp;
      if (re_tmp == 0.0) {
        u0 = k - hnRows;
        ytmp[u0].re = ar / static_cast<double>(hnRows);
        ytmp[u0].im = 0.0;
      } else if (ar == 0.0) {
        u0 = k - hnRows;
        ytmp[u0].re = 0.0;
        ytmp[u0].im = re_tmp / static_cast<double>(hnRows);
      } else {
        u0 = k - hnRows;
        ytmp[u0].re = ar / static_cast<double>(hnRows);
        ytmp[u0].im = re_tmp / static_cast<double>(hnRows);
      }
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        dc, ar_tmp, b_ar_tmp, c_ar_tmp, d_ar_tmp, b_ar)

    for (int k = hnRows; k <= minHnrowsNxBy2; k++) {
      ar_tmp = wwc[k - 1].re;
      b_ar_tmp = fv[k - 1].im;
      c_ar_tmp = wwc[k - 1].im;
      d_ar_tmp = fv[k - 1].re;
      b_ar = ar_tmp * d_ar_tmp + c_ar_tmp * b_ar_tmp;
      ar_tmp = ar_tmp * b_ar_tmp - c_ar_tmp * d_ar_tmp;
      if (ar_tmp == 0.0) {
        dc.re = b_ar / static_cast<double>(hnRows);
        dc.im = 0.0;
      } else if (b_ar == 0.0) {
        dc.re = 0.0;
        dc.im = ar_tmp / static_cast<double>(hnRows);
      } else {
        dc.re = b_ar / static_cast<double>(hnRows);
        dc.im = ar_tmp / static_cast<double>(hnRows);
      }
      ytmp[k - hnRows] = dc;
    }
  }
  if (static_cast<int>(hnRows < 1600)) {
    for (int c_i{0}; c_i < hnRows; c_i++) {
      double b_ytmp_re_tmp;
      double ytmp_im;
      double ytmp_re_tmp;
      i2 = wrapIndex[c_i];
      re_tmp = ytmp[c_i].re;
      b_re_tmp = reconVar1[c_i].im;
      c_re_tmp = ytmp[c_i].im;
      e_re_tmp = reconVar1[c_i].re;
      ar = ytmp[i2 - 1].re;
      ytmp_im = -ytmp[i2 - 1].im;
      ytmp_re_tmp = reconVar2[c_i].im;
      b_ytmp_re_tmp = reconVar2[c_i].re;
      y[c_i].re = 0.5 * ((re_tmp * e_re_tmp - c_re_tmp * b_re_tmp) +
                         (ar * b_ytmp_re_tmp - ytmp_im * ytmp_re_tmp));
      y[c_i].im = 0.5 * ((re_tmp * b_re_tmp + c_re_tmp * e_re_tmp) +
                         (ar * ytmp_re_tmp + ytmp_im * b_ytmp_re_tmp));
      u0 = hnRows + c_i;
      y[u0].re = 0.5 * ((re_tmp * b_ytmp_re_tmp - c_re_tmp * ytmp_re_tmp) +
                        (ar * e_re_tmp - ytmp_im * b_re_tmp));
      y[u0].im = 0.5 * ((re_tmp * ytmp_re_tmp + c_re_tmp * b_ytmp_re_tmp) +
                        (ar * b_re_tmp + ytmp_im * e_re_tmp));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        i2, c_ytmp_re_tmp, d_ytmp_re_tmp, e_ytmp_re_tmp, f_ytmp_re_tmp,        \
            ytmp_re, b_ytmp_im, g_ytmp_re_tmp, h_ytmp_re_tmp)

    for (int c_i = 0; c_i < hnRows; c_i++) {
      i2 = wrapIndex[c_i];
      c_ytmp_re_tmp = ytmp[c_i].re;
      d_ytmp_re_tmp = reconVar1[c_i].im;
      e_ytmp_re_tmp = ytmp[c_i].im;
      f_ytmp_re_tmp = reconVar1[c_i].re;
      ytmp_re = ytmp[i2 - 1].re;
      b_ytmp_im = -ytmp[i2 - 1].im;
      g_ytmp_re_tmp = reconVar2[c_i].im;
      h_ytmp_re_tmp = reconVar2[c_i].re;
      y[c_i].re =
          0.5 *
          ((c_ytmp_re_tmp * f_ytmp_re_tmp - e_ytmp_re_tmp * d_ytmp_re_tmp) +
           (ytmp_re * h_ytmp_re_tmp - b_ytmp_im * g_ytmp_re_tmp));
      y[c_i].im =
          0.5 *
          ((c_ytmp_re_tmp * d_ytmp_re_tmp + e_ytmp_re_tmp * f_ytmp_re_tmp) +
           (ytmp_re * g_ytmp_re_tmp + b_ytmp_im * h_ytmp_re_tmp));
      i2 = hnRows + c_i;
      y[i2].re =
          0.5 *
          ((c_ytmp_re_tmp * h_ytmp_re_tmp - e_ytmp_re_tmp * g_ytmp_re_tmp) +
           (ytmp_re * f_ytmp_re_tmp - b_ytmp_im * d_ytmp_re_tmp));
      y[i2].im =
          0.5 *
          ((c_ytmp_re_tmp * g_ytmp_re_tmp + e_ytmp_re_tmp * h_ytmp_re_tmp) +
           (ytmp_re * d_ytmp_re_tmp + b_ytmp_im * f_ytmp_re_tmp));
    }
  }
  u0 = y.size(0);
  if (static_cast<int>(u0 < 1600)) {
    for (int i3{0}; i3 < u0; i3++) {
      re_tmp = y[i3].re;
      c_re_tmp = y[i3].im;
      if (c_re_tmp == 0.0) {
        b_re_tmp = re_tmp / 2.0;
        re_tmp = 0.0;
      } else if (re_tmp == 0.0) {
        b_re_tmp = 0.0;
        re_tmp = c_re_tmp / 2.0;
      } else {
        b_re_tmp = re_tmp / 2.0;
        re_tmp = c_re_tmp / 2.0;
      }
      y[i3].re = b_re_tmp;
      y[i3].im = re_tmp;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(im, ai, re)

    for (int i3 = 0; i3 < u0; i3++) {
      im = y[i3].re;
      ai = y[i3].im;
      if (ai == 0.0) {
        re = im / 2.0;
        im = 0.0;
      } else if (im == 0.0) {
        re = 0.0;
        im = ai / 2.0;
      } else {
        re = im / 2.0;
        im = ai / 2.0;
      }
      y[i3].re = re;
      y[i3].im = im;
    }
  }
}

//
// Arguments    : const array<double, 1U> &x
//                array<creal_T, 1U> &y
//                int nrowsx
//                int nRows
//                int nfft
//                const array<creal_T, 1U> &wwc
//                const array<double, 2U> &costab
//                const array<double, 2U> &sintab
//                const array<double, 2U> &costabinv
//                const array<double, 2U> &sintabinv
// Return Type  : void
//
void FFTImplementationCallback::doHalfLengthBluestein(
    const array<double, 1U> &x, array<creal_T, 1U> &y, int nrowsx, int nRows,
    int nfft, const array<creal_T, 1U> &wwc, const array<double, 2U> &costab,
    const array<double, 2U> &sintab, const array<double, 2U> &costabinv,
    const array<double, 2U> &sintabinv)
{
  array<creal_T, 1U> b_fv;
  array<creal_T, 1U> fv;
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<creal_T, 1U> ytmp;
  array<double, 2U> a__1;
  array<double, 2U> costable;
  array<double, 2U> hcostabinv;
  array<double, 2U> hsintab;
  array<double, 2U> hsintabinv;
  array<double, 2U> sintable;
  array<int, 2U> wrapIndex;
  creal_T a;
  creal_T b;
  double b_re_tmp;
  double b_ytmp_im;
  double b_ytmp_re;
  double c_re_tmp;
  double c_ytmp_re_tmp;
  double d;
  double d1;
  double d2;
  double d3;
  double d_re_tmp;
  double d_ytmp_re_tmp;
  double e_re_tmp;
  double e_ytmp_re_tmp;
  double f_re_tmp;
  double f_ytmp_re_tmp;
  double g_re_tmp;
  double g_ytmp_re_tmp;
  double h_re_tmp;
  double h_ytmp_re_tmp;
  double re_tmp;
  int b_tmp;
  int hnRows;
  int i2;
  int i3;
  int minHnrowsNxBy2;
  int u0;
  boolean_T nxeven;
  hnRows = nRows / 2;
  ytmp.set_size(hnRows);
  if (hnRows > nrowsx) {
    ytmp.set_size(hnRows);
    if (hnRows - 1 >= 0) {
      std::memset(&ytmp[0], 0,
                  static_cast<unsigned int>(hnRows) * sizeof(creal_T));
    }
  }
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    nxeven = true;
    u0 = x.size(0);
  } else if (x.size(0) >= nRows) {
    nxeven = true;
    u0 = nRows;
  } else {
    nxeven = false;
    u0 = x.size(0) - 1;
  }
  FFTImplementationCallback::generate_twiddle_tables(nRows << 1, costable,
                                                     sintable, a__1);
  FFTImplementationCallback::get_half_twiddle_tables(costab, sintab, costabinv,
                                                     sintabinv, a__1, hsintab,
                                                     hcostabinv, hsintabinv);
  reconVar1.set_size(hnRows);
  reconVar2.set_size(hnRows);
  wrapIndex.set_size(1, hnRows);
  for (int i{0}; i < hnRows; i++) {
    minHnrowsNxBy2 = i << 1;
    re_tmp = sintable[minHnrowsNxBy2];
    b_re_tmp = costable[minHnrowsNxBy2];
    reconVar1[i].re = re_tmp + 1.0;
    reconVar1[i].im = -b_re_tmp;
    reconVar2[i].re = 1.0 - re_tmp;
    reconVar2[i].im = b_re_tmp;
    if (i != 0) {
      wrapIndex[i] = (hnRows - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
  }
  if (u0 > nRows) {
    u0 = nRows;
  }
  minHnrowsNxBy2 = u0 / 2 - 1;
  if (static_cast<int>(minHnrowsNxBy2 + 1 < 1600)) {
    for (int k1{0}; k1 <= minHnrowsNxBy2; k1++) {
      a = wwc[(hnRows + k1) - 1];
      u0 = k1 << 1;
      b.re = x[u0];
      b.im = x[u0 + 1];
      ytmp[k1].re = a.re * b.re + a.im * b.im;
      ytmp[k1].im = a.re * b.im - a.im * b.re;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b, a, b_tmp)

    for (int k1 = 0; k1 <= minHnrowsNxBy2; k1++) {
      a = wwc[(hnRows + k1) - 1];
      b_tmp = k1 << 1;
      b.re = x[b_tmp];
      b.im = x[b_tmp + 1];
      ytmp[k1].re = a.re * b.re + a.im * b.im;
      ytmp[k1].im = a.re * b.im - a.im * b.re;
    }
  }
  if (!nxeven) {
    a = wwc[hnRows + minHnrowsNxBy2];
    if (minHnrowsNxBy2 < 0) {
      u0 = 0;
    } else {
      u0 = (minHnrowsNxBy2 + 1) << 1;
    }
    b.re = x[u0];
    ytmp[minHnrowsNxBy2 + 1].re = a.re * b.re + a.im * 0.0;
    ytmp[minHnrowsNxBy2 + 1].im = a.re * 0.0 - a.im * b.re;
    if (minHnrowsNxBy2 + 3 <= hnRows) {
      u0 = minHnrowsNxBy2 + 3;
      if (u0 <= hnRows) {
        std::memset(&ytmp[u0 + -1], 0,
                    static_cast<unsigned int>((hnRows - u0) + 1) *
                        sizeof(creal_T));
      }
    }
  } else if (minHnrowsNxBy2 + 2 <= hnRows) {
    u0 = minHnrowsNxBy2 + 2;
    if (u0 <= hnRows) {
      std::memset(&ytmp[u0 + -1], 0,
                  static_cast<unsigned int>((hnRows - u0) + 1) *
                      sizeof(creal_T));
    }
  }
  minHnrowsNxBy2 = nfft / 2;
  FFTImplementationCallback::r2br_r2dit_trig_impl(ytmp, minHnrowsNxBy2, a__1,
                                                  hsintab, fv);
  FFTImplementationCallback::r2br_r2dit_trig_impl(wwc, minHnrowsNxBy2, a__1,
                                                  hsintab, b_fv);
  u0 = fv.size(0);
  b_fv.set_size(fv.size(0));
  if (static_cast<int>(fv.size(0) < 1600)) {
    for (int b_i{0}; b_i < u0; b_i++) {
      re_tmp = fv[b_i].re;
      b_re_tmp = b_fv[b_i].im;
      c_re_tmp = fv[b_i].im;
      e_re_tmp = b_fv[b_i].re;
      b_fv[b_i].re = re_tmp * e_re_tmp - c_re_tmp * b_re_tmp;
      b_fv[b_i].im = re_tmp * b_re_tmp + c_re_tmp * e_re_tmp;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d_re_tmp, f_re_tmp, g_re_tmp, h_re_tmp)

    for (int b_i = 0; b_i < u0; b_i++) {
      d_re_tmp = fv[b_i].re;
      f_re_tmp = b_fv[b_i].im;
      g_re_tmp = fv[b_i].im;
      h_re_tmp = b_fv[b_i].re;
      b_fv[b_i].re = d_re_tmp * h_re_tmp - g_re_tmp * f_re_tmp;
      b_fv[b_i].im = d_re_tmp * f_re_tmp + g_re_tmp * h_re_tmp;
    }
  }
  FFTImplementationCallback::r2br_r2dit_trig_impl(b_fv, minHnrowsNxBy2,
                                                  hcostabinv, hsintabinv, fv);
  if (fv.size(0) > 1) {
    re_tmp = 1.0 / static_cast<double>(fv.size(0));
    u0 = fv.size(0);
    if (static_cast<int>(fv.size(0) < 1600)) {
      for (int i1{0}; i1 < u0; i1++) {
        fv[i1].re = re_tmp * fv[i1].re;
        fv[i1].im = re_tmp * fv[i1].im;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (int i1 = 0; i1 < u0; i1++) {
        fv[i1].re = re_tmp * fv[i1].re;
        fv[i1].im = re_tmp * fv[i1].im;
      }
    }
  }
  minHnrowsNxBy2 = wwc.size(0);
  if (static_cast<int>((wwc.size(0) - hnRows) + 1 < 1600)) {
    for (int k{hnRows}; k <= minHnrowsNxBy2; k++) {
      re_tmp = wwc[k - 1].re;
      b_re_tmp = fv[k - 1].im;
      c_re_tmp = wwc[k - 1].im;
      e_re_tmp = fv[k - 1].re;
      u0 = k - hnRows;
      ytmp[u0].re = re_tmp * e_re_tmp + c_re_tmp * b_re_tmp;
      ytmp[u0].im = re_tmp * b_re_tmp - c_re_tmp * e_re_tmp;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d, d1, d2, d3, i3)

    for (int k = hnRows; k <= minHnrowsNxBy2; k++) {
      d = wwc[k - 1].re;
      d1 = fv[k - 1].im;
      d2 = wwc[k - 1].im;
      d3 = fv[k - 1].re;
      i3 = k - hnRows;
      ytmp[i3].re = d * d3 + d2 * d1;
      ytmp[i3].im = d * d1 - d2 * d3;
    }
  }
  if (static_cast<int>(hnRows < 1600)) {
    for (int c_i{0}; c_i < hnRows; c_i++) {
      double b_ytmp_re_tmp;
      double ytmp_im;
      double ytmp_re;
      double ytmp_re_tmp;
      i2 = wrapIndex[c_i];
      re_tmp = ytmp[c_i].re;
      b_re_tmp = reconVar1[c_i].im;
      c_re_tmp = ytmp[c_i].im;
      e_re_tmp = reconVar1[c_i].re;
      ytmp_re = ytmp[i2 - 1].re;
      ytmp_im = -ytmp[i2 - 1].im;
      ytmp_re_tmp = reconVar2[c_i].im;
      b_ytmp_re_tmp = reconVar2[c_i].re;
      y[c_i].re = 0.5 * ((re_tmp * e_re_tmp - c_re_tmp * b_re_tmp) +
                         (ytmp_re * b_ytmp_re_tmp - ytmp_im * ytmp_re_tmp));
      y[c_i].im = 0.5 * ((re_tmp * b_re_tmp + c_re_tmp * e_re_tmp) +
                         (ytmp_re * ytmp_re_tmp + ytmp_im * b_ytmp_re_tmp));
      u0 = hnRows + c_i;
      y[u0].re = 0.5 * ((re_tmp * b_ytmp_re_tmp - c_re_tmp * ytmp_re_tmp) +
                        (ytmp_re * e_re_tmp - ytmp_im * b_re_tmp));
      y[u0].im = 0.5 * ((re_tmp * ytmp_re_tmp + c_re_tmp * b_ytmp_re_tmp) +
                        (ytmp_re * b_re_tmp + ytmp_im * e_re_tmp));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        i2, c_ytmp_re_tmp, d_ytmp_re_tmp, e_ytmp_re_tmp, f_ytmp_re_tmp,        \
            b_ytmp_re, b_ytmp_im, g_ytmp_re_tmp, h_ytmp_re_tmp)

    for (int c_i = 0; c_i < hnRows; c_i++) {
      i2 = wrapIndex[c_i];
      c_ytmp_re_tmp = ytmp[c_i].re;
      d_ytmp_re_tmp = reconVar1[c_i].im;
      e_ytmp_re_tmp = ytmp[c_i].im;
      f_ytmp_re_tmp = reconVar1[c_i].re;
      b_ytmp_re = ytmp[i2 - 1].re;
      b_ytmp_im = -ytmp[i2 - 1].im;
      g_ytmp_re_tmp = reconVar2[c_i].im;
      h_ytmp_re_tmp = reconVar2[c_i].re;
      y[c_i].re =
          0.5 *
          ((c_ytmp_re_tmp * f_ytmp_re_tmp - e_ytmp_re_tmp * d_ytmp_re_tmp) +
           (b_ytmp_re * h_ytmp_re_tmp - b_ytmp_im * g_ytmp_re_tmp));
      y[c_i].im =
          0.5 *
          ((c_ytmp_re_tmp * d_ytmp_re_tmp + e_ytmp_re_tmp * f_ytmp_re_tmp) +
           (b_ytmp_re * g_ytmp_re_tmp + b_ytmp_im * h_ytmp_re_tmp));
      i2 = hnRows + c_i;
      y[i2].re =
          0.5 *
          ((c_ytmp_re_tmp * h_ytmp_re_tmp - e_ytmp_re_tmp * g_ytmp_re_tmp) +
           (b_ytmp_re * f_ytmp_re_tmp - b_ytmp_im * d_ytmp_re_tmp));
      y[i2].im =
          0.5 *
          ((c_ytmp_re_tmp * g_ytmp_re_tmp + e_ytmp_re_tmp * h_ytmp_re_tmp) +
           (b_ytmp_re * d_ytmp_re_tmp + b_ytmp_im * f_ytmp_re_tmp));
    }
  }
}

//
// Arguments    : int nfft
//                boolean_T useRadix2
//                int &nRows
// Return Type  : int
//
int FFTImplementationCallback::get_algo_sizes(int nfft, boolean_T useRadix2,
                                              int &nRows)
{
  int n2blue;
  n2blue = 1;
  if (useRadix2) {
    nRows = nfft;
  } else {
    if (nfft > 0) {
      int pmax;
      n2blue = (nfft + nfft) - 1;
      pmax = 31;
      if (n2blue <= 1) {
        pmax = 0;
      } else {
        int pmin;
        boolean_T exitg1;
        pmin = 0;
        exitg1 = false;
        while ((!exitg1) && (pmax - pmin > 1)) {
          int k;
          int pow2p;
          k = (pmin + pmax) >> 1;
          pow2p = 1 << k;
          if (pow2p == n2blue) {
            pmax = k;
            exitg1 = true;
          } else if (pow2p > n2blue) {
            pmax = k;
          } else {
            pmin = k;
          }
        }
      }
      n2blue = 1 << pmax;
    }
    nRows = n2blue;
  }
  return n2blue;
}

//
// Arguments    : const array<creal_T, 1U> &x
//                int unsigned_nRows
//                const array<double, 2U> &costab
//                const array<double, 2U> &sintab
//                array<creal_T, 1U> &y
// Return Type  : void
//
void FFTImplementationCallback::r2br_r2dit_trig_impl(
    const array<creal_T, 1U> &x, int unsigned_nRows,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    array<creal_T, 1U> &y)
{
  double im;
  double re;
  double temp_im;
  double temp_re;
  int iDelta;
  int iDelta2;
  int iheight;
  int istart;
  int iy;
  int j;
  int ju;
  int k;
  int nRowsD2;
  y.set_size(unsigned_nRows);
  if (unsigned_nRows > x.size(0)) {
    y.set_size(unsigned_nRows);
    for (int i{0}; i < unsigned_nRows; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  j = x.size(0);
  if (j > unsigned_nRows) {
    j = unsigned_nRows;
  }
  iDelta = unsigned_nRows - 2;
  nRowsD2 = unsigned_nRows / 2;
  k = nRowsD2 / 2;
  iy = 0;
  ju = 0;
  for (int i{0}; i <= j - 2; i++) {
    boolean_T tst;
    y[iy] = x[i];
    istart = unsigned_nRows;
    tst = true;
    while (tst) {
      istart >>= 1;
      ju ^= istart;
      tst = ((ju & istart) == 0);
    }
    iy = ju;
  }
  if (j - 2 < 0) {
    istart = 0;
  } else {
    istart = j - 1;
  }
  y[iy] = x[istart];
  if (unsigned_nRows > 1) {
    for (int i{0}; i <= iDelta; i += 2) {
      temp_re = y[i + 1].re;
      temp_im = y[i + 1].im;
      re = y[i].re;
      im = y[i].im;
      y[i + 1].re = re - temp_re;
      y[i + 1].im = y[i].im - y[i + 1].im;
      re += temp_re;
      im += temp_im;
      y[i].re = re;
      y[i].im = im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_i;
    for (b_i = 0; b_i < iheight; b_i += iDelta2) {
      istart = b_i + iDelta;
      temp_re = y[istart].re;
      temp_im = y[istart].im;
      y[istart].re = y[b_i].re - temp_re;
      y[istart].im = y[b_i].im - temp_im;
      y[b_i].re = y[b_i].re + temp_re;
      y[b_i].im = y[b_i].im + temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      double twid_im;
      double twid_re;
      twid_re = costab[j];
      twid_im = sintab[j];
      b_i = istart;
      iy = istart + iheight;
      while (b_i < iy) {
        ju = b_i + iDelta;
        re = y[ju].im;
        im = y[ju].re;
        temp_re = twid_re * im - twid_im * re;
        temp_im = twid_re * re + twid_im * im;
        y[ju].re = y[b_i].re - temp_re;
        y[ju].im = y[b_i].im - temp_im;
        y[b_i].re = y[b_i].re + temp_re;
        y[b_i].im = y[b_i].im + temp_im;
        b_i += iDelta2;
      }
      istart++;
    }
    k = static_cast<int>(static_cast<unsigned int>(k) >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
}

//
// Arguments    : const array<double, 1U> &x
//                int unsigned_nRows
//                const array<double, 2U> &costab
//                const array<double, 2U> &sintab
//                array<creal_T, 1U> &y
// Return Type  : void
//
void FFTImplementationCallback::r2br_r2dit_trig_impl(
    const array<double, 1U> &x, int unsigned_nRows,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    array<creal_T, 1U> &y)
{
  y.set_size(unsigned_nRows);
  if (unsigned_nRows > x.size(0)) {
    y.set_size(unsigned_nRows);
    for (int i{0}; i < unsigned_nRows; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  if (unsigned_nRows != 1) {
    FFTImplementationCallback::doHalfLengthRadix2(x, y, unsigned_nRows, costab,
                                                  sintab);
  } else {
    y[0].re = x[0];
    y[0].im = 0.0;
  }
}

} // namespace fft
} // namespace internal
} // namespace coder

//
// File trailer for FFTImplementationCallback.cpp
//
// [EOF]
//
