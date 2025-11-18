//
// File: swsmodel.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "swsmodel.h"
#include "abs.h"
#include "angle.h"
#include "anynan.h"
#include "filter.h"
#include "hanning.h"
#include "inv.h"
#include "length.h"
#include "lpcfit.h"
#include "mean.h"
#include "padarray.h"
#include "resample.h"
#include "roots.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "swsmodel_data.h"
#include "swsmodel_initialize.h"
#include "toeplitz.h"
#include "xcorr.h"
#include "coder_array.h"
#include "omp.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <emmintrin.h>

// Function Definitions
//
// [F,M] = swsmodel(D,R,H)  Sine wave speech analysis
//        D is a speech example sampled at R samples per second.
//        Return a sinusoid model of up to 4 components with each sinusoid
//        defined by a row of F (frequencies in Hz) and M (linear magnitude).
//        Each column of F and M corresponds to H samples at R.
//        Rows of F are sorted with lowest frequency first;
//        sinusoids cannot cross.
//
//        Relies on lpcfit.m and lpca2frq.m to form LPC model and convert it
//        into frequencies.
//  2001-03-12 dpwe@ee.columbia.edu $Header: $
//
// Arguments    : coder::array<double, 2U> &D
//                double R
//                double H
//                coder::array<double, 2U> &F
//                coder::array<double, 2U> &M
// Return Type  : void
//
void swsmodel(coder::array<double, 2U> &D, double R, double H,
              coder::array<double, 2U> &F, coder::array<double, 2U> &M)
{
  __m128d r1;
  coder::array<double, 2U> b_D;
  coder::array<double, 2U> fa;
  coder::array<double, 2U> lpca;
  coder::array<double, 2U> ma;
  coder::array<double, 2U> rs;
  coder::array<double, 2U> wxx;
  coder::array<double, 2U> xx;
  coder::array<double, 1U> g;
  coder::array<double, 1U> r;
  creal_T rts_data[32];
  creal_T b_rts_data[8];
  double dv2[64];
  double frqs_data[32];
  double mags_data[32];
  double y_data[32];
  double rxx[9];
  double b_lpca[8];
  double dv[2];
  double HH;
  double b_hop;
  double nhops;
  double w;
  int iidx_data[32];
  int frqs_size[2];
  int rts_size[2];
  int b_loop_ub;
  int c_loop_ub;
  int d_loop_ub;
  int i1;
  int loop_ub;
  int vectorUB;
  signed char tmp_data[32];
  if (!isInitialized_swsmodel) {
    swsmodel_initialize();
  }
  //  Target sampling rate
  //  Resample to 8 kHz, so LPC only picks main formants
  if (R != 6000.0) {
    b_D.set_size(D.size(0), D.size(1));
    loop_ub = D.size(0) * D.size(1) - 1;
    for (int i{0}; i <= loop_ub; i++) {
      b_D[i] = D[i];
    }
    coder::resample(b_D, std::round(R / 1000.0), D);
  }
  //  Step size in units of my sampling rate (force to be even integer)
  HH = 2.0 * std::round(H / R * 6000.0 / 2.0);
  //  Form 8th-order LPC model (3 or 4 pole pairs)
  w = 2.0 * HH;
  //  [a,g,e] = lpcfit(x,p,h,w,ov)  Fit LPC to short-time segments
  //     x is a stretch of signal.  Using w point (2*h) windows every
  //     h points (128), fit order p LPC models.  Return the successive
  //     all-pole coefficients as rows of a, the per-frame gains in g
  //     and the residual excitation in e.
  //     ov nonzero selects overlap-add of window-length
  //     residuals, otherwise successive hop-sized residuals are concatenated
  //     for independent near-perfect reconstruction with lpcsynth.
  //     (default is 1)
  //  2001-02-25 dpwe@ee.columbia.edu $Header:
  //  /homes/dpwe/matlab/columbiafns/RCS/lpcfit.m,v 1.1 2004/03/30 20:55:52 dpwe
  //  Exp $
  if (D.size(1) == 1) {
    loop_ub = D.size(0);
    b_D.set_size(1, D.size(0));
    for (int i{0}; i < loop_ub; i++) {
      b_D[b_D.size(0) * i] = D[i];
    }
    D.set_size(1, loop_ub);
    for (int i{0}; i < loop_ub; i++) {
      D[i] = b_D[i];
    }
    //  Convert X from column to row
  }
  nhops = std::floor(
      static_cast<double>(coder::internal::intlength(D.size(0), D.size(1))) /
      HH);
  //  Pad x with zeros so that we can extract complete w-length windows
  //  from it
  b_D.set_size(D.size(0), D.size(1));
  loop_ub = D.size(0) * D.size(1) - 1;
  for (int i{0}; i <= loop_ub; i++) {
    b_D[i] = D[i];
  }
  coder::padarray(b_D, (w - HH) / 2.0, D);
  b_loop_ub = static_cast<int>(nhops);
  lpca.set_size(static_cast<int>(nhops), 9);
  g.set_size(static_cast<int>(nhops));
  //  Pre-emphasis
  b_D.set_size(D.size(0), D.size(1));
  loop_ub = D.size(0) * D.size(1) - 1;
  for (int i{0}; i <= loop_ub; i++) {
    b_D[i] = D[i];
  }
  coder::filter(b_D, D);
  for (int hop{0}; hop < b_loop_ub; hop++) {
    __m128d r2;
    //  Extract segment of signal
    if (std::isnan(w)) {
      wxx.set_size(1, 1);
      wxx[0] = rtNaN;
    } else if (w < 1.0) {
      wxx.set_size(1, 0);
    } else {
      wxx.set_size(1, static_cast<int>(w - 1.0) + 1);
      loop_ub = static_cast<int>(w - 1.0);
      vectorUB = ((static_cast<int>(w - 1.0) + 1) / 2) << 1;
      c_loop_ub = vectorUB - 2;
      for (int i{0}; i <= c_loop_ub; i += 2) {
        dv[0] = i;
        dv[1] = i + 1;
        r1 = _mm_loadu_pd(&dv[0]);
        _mm_storeu_pd(&wxx[i], _mm_add_pd(_mm_set1_pd(1.0), r1));
      }
      for (int i{vectorUB}; i <= loop_ub; i++) {
        wxx[i] = static_cast<double>(i) + 1.0;
      }
    }
    c_loop_ub = wxx.size(1);
    xx.set_size(1, wxx.size(1));
    b_hop = ((static_cast<double>(hop) + 1.0) - 1.0) * HH;
    for (int i{0}; i < c_loop_ub; i++) {
      xx[i] = D[static_cast<int>(b_hop + wxx[i]) - 1];
    }
    //  Apply hanning window
    coder::hanning(w, r);
    if (r.size(0) == xx.size(1)) {
      wxx.set_size(1, wxx.size(1));
      loop_ub = (xx.size(1) / 2) << 1;
      vectorUB = loop_ub - 2;
      for (int i{0}; i <= vectorUB; i += 2) {
        r1 = _mm_loadu_pd(&xx[i]);
        r2 = _mm_loadu_pd(&r[i]);
        _mm_storeu_pd(&wxx[i], _mm_mul_pd(r1, r2));
      }
      for (int i{loop_ub}; i < c_loop_ub; i++) {
        wxx[i] = xx[i] * r[i];
      }
    } else {
      binary_expand_op(wxx, xx, r);
    }
    //  Form autocorrelation (calculates *way* too many points)
    coder::xcorr(wxx, xx);
    //  extract just the points we need (middle p+1 points)
    for (int i{0}; i < 9; i++) {
      rxx[i] = xx[static_cast<int>(w + static_cast<double>(i)) - 1];
    }
    double dv1[64];
    //  Setup the normal equations
    //  Solve for a (horribly inefficient to use full inv())
    //  Calculate residual by filtering windowed xx
    coder::toeplitz(&rxx[0], dv1);
    coder::inv(dv1, dv2);
    std::memset(&b_lpca[0], 0, sizeof(double) << 3);
    for (int i{0}; i < 8; i++) {
      __m128d r3;
      loop_ub = i << 3;
      r1 = _mm_loadu_pd(&dv2[loop_ub]);
      r2 = _mm_loadu_pd(&b_lpca[0]);
      r3 = _mm_set1_pd(rxx[i + 1]);
      _mm_storeu_pd(&b_lpca[0], _mm_add_pd(r2, _mm_mul_pd(r1, r3)));
      r1 = _mm_loadu_pd(&dv2[loop_ub + 2]);
      r2 = _mm_loadu_pd(&b_lpca[2]);
      _mm_storeu_pd(&b_lpca[2], _mm_add_pd(r2, _mm_mul_pd(r1, r3)));
      r1 = _mm_loadu_pd(&dv2[loop_ub + 4]);
      r2 = _mm_loadu_pd(&b_lpca[4]);
      _mm_storeu_pd(&b_lpca[4], _mm_add_pd(r2, _mm_mul_pd(r1, r3)));
      r1 = _mm_loadu_pd(&dv2[loop_ub + 6]);
      r2 = _mm_loadu_pd(&b_lpca[6]);
      _mm_storeu_pd(&b_lpca[6], _mm_add_pd(r2, _mm_mul_pd(r1, r3)));
    }
    rxx[0] = 1.0;
    r1 = _mm_loadu_pd(&b_lpca[0]);
    r2 = _mm_set1_pd(-1.0);
    _mm_storeu_pd(&rxx[1], _mm_mul_pd(r1, r2));
    r1 = _mm_loadu_pd(&b_lpca[2]);
    _mm_storeu_pd(&rxx[3], _mm_mul_pd(r1, r2));
    r1 = _mm_loadu_pd(&b_lpca[4]);
    _mm_storeu_pd(&rxx[5], _mm_mul_pd(r1, r2));
    r1 = _mm_loadu_pd(&b_lpca[6]);
    _mm_storeu_pd(&rxx[7], _mm_mul_pd(r1, r2));
    coder::filter(rxx, wxx, rs);
    loop_ub = rs.size(1);
    xx.set_size(1, rs.size(1));
    for (int i{0}; i < loop_ub; i++) {
      b_hop = rs[i];
      xx[i] = b_hop * b_hop;
    }
    //  Save filter, gain and residual
    for (int i{0}; i < 9; i++) {
      lpca[hop + lpca.size(0) * i] = rxx[i];
    }
    g[hop] = std::sqrt(coder::mean(xx));
  }
  //  Throw away first (win-hop)/2 pts if in overlap mode
  //  for proper synchronization of resynth
  //  Convert poles to sorted freqs and magnitudes
  //  If only 3 nonzero freqs are found, 4th row will have mag/frq zero
  //  [f,m] = lpca2frq(a,g)  Convert LPC analysis frames into resonant
  //  frequencies
  //     Each row of a defines an LPC analysis e.g. from lpcfit.  Convert
  //     this into poles, and return the frequencies in rows of f.
  //  2001-02-25 dpwe@ee.columbia.edu
  fa.set_size(static_cast<int>(nhops), 4);
  d_loop_ub = lpca.size(0) << 2;
  if (d_loop_ub - 1 >= 0) {
    std::memset(&fa[0], 0,
                static_cast<unsigned int>(d_loop_ub) * sizeof(double));
  }
  ma.set_size(static_cast<int>(nhops), 4);
  if (d_loop_ub - 1 >= 0) {
    std::memset(&ma[0], 0,
                static_cast<unsigned int>(d_loop_ub) * sizeof(double));
  }
  for (int c_hop{0}; c_hop < b_loop_ub; c_hop++) {
    creal_T c_rts_data[32];
    int rts_size_idx_1;
    for (int i{0}; i < 8; i++) {
      b_lpca[i] = lpca[c_hop + lpca.size(0) * (i + 1)];
    }
    if (coder::anynan(b_lpca)) {
      c_loop_ub = 4;
      rts_size_idx_1 = 4;
      std::memset(&rts_data[0], 0, 16U * sizeof(creal_T));
    } else {
      rxx[0] = 1.0;
      for (int i{0}; i < 8; i++) {
        rxx[i + 1] = lpca[c_hop + lpca.size(0) * (i + 1)];
      }
      c_loop_ub = coder::roots(rxx, b_rts_data);
      rts_size_idx_1 = 1;
      if (c_loop_ub - 1 >= 0) {
        std::copy(&b_rts_data[0], &b_rts_data[c_loop_ub], &rts_data[0]);
      }
    }
    rts_size[0] = rts_size_idx_1;
    rts_size[1] = c_loop_ub;
    for (int i{0}; i < c_loop_ub; i++) {
      for (int hop{0}; hop < rts_size_idx_1; hop++) {
        loop_ub = i + c_loop_ub * hop;
        vectorUB = hop + rts_size_idx_1 * i;
        c_rts_data[vectorUB].re = rts_data[loop_ub].re;
        c_rts_data[vectorUB].im = -rts_data[loop_ub].im;
      }
    }
    coder::angle(c_rts_data, rts_size, frqs_data, frqs_size);
    for (int i{0}; i < c_loop_ub; i++) {
      for (int hop{0}; hop < rts_size_idx_1; hop++) {
        loop_ub = i + c_loop_ub * hop;
        vectorUB = hop + rts_size_idx_1 * i;
        c_rts_data[vectorUB].re = rts_data[loop_ub].re;
        c_rts_data[vectorUB].im = -rts_data[loop_ub].im;
      }
    }
    int mags_size[2];
    coder::b_abs(c_rts_data, rts_size, mags_data, mags_size);
    loop_ub = mags_size[0] * mags_size[1];
    b_hop = g[c_hop];
    vectorUB = (loop_ub / 2) << 1;
    c_loop_ub = vectorUB - 2;
    for (int i{0}; i <= c_loop_ub; i += 2) {
      r1 = _mm_loadu_pd(&mags_data[i]);
      _mm_storeu_pd(
          &mags_data[i],
          _mm_div_pd(_mm_set1_pd(b_hop), _mm_sub_pd(_mm_set1_pd(1.0), r1)));
    }
    for (int i{vectorUB}; i < loop_ub; i++) {
      mags_data[i] = b_hop / (1.0 - mags_data[i]);
    }
    loop_ub = frqs_size[0] * frqs_size[1];
    if (loop_ub - 1 >= 0) {
      std::copy(&frqs_data[0], &frqs_data[loop_ub], &y_data[0]);
    }
    coder::internal::sort(y_data, frqs_size, iidx_data, rts_size);
    loop_ub = rts_size[0] * rts_size[1];
    rts_size_idx_1 = 0;
    vectorUB = 0;
    for (int i{0}; i < loop_ub; i++) {
      c_loop_ub = iidx_data[i];
      y_data[i] = c_loop_ub;
      if (frqs_data[c_loop_ub - 1] > 0.0) {
        rts_size_idx_1++;
        tmp_data[vectorUB] = static_cast<signed char>(i);
        vectorUB++;
      }
    }
    for (int i{0}; i < rts_size_idx_1; i++) {
      loop_ub = static_cast<signed char>(y_data[tmp_data[i]]) - 1;
      fa[c_hop + fa.size(0) * i] = frqs_data[loop_ub];
      ma[c_hop + ma.size(0) * i] = mags_data[loop_ub];
    }
  }
  //  Convert frqs into Hz
  F.set_size(4, static_cast<int>(nhops));
  loop_ub = static_cast<int>(nhops);
  M.set_size(4, static_cast<int>(nhops));
  if (static_cast<int>(d_loop_ub < 1600)) {
    for (int b_i{0}; b_i < b_loop_ub; b_i++) {
      F[4 * b_i] = fa[b_i] * 6000.0 / 6.2831853071795862;
      M[4 * b_i] = ma[b_i];
      loop_ub = 4 * b_i + 1;
      F[loop_ub] = fa[b_i + fa.size(0)] * 6000.0 / 6.2831853071795862;
      M[loop_ub] = ma[b_i + ma.size(0)];
      loop_ub = 4 * b_i + 2;
      F[loop_ub] = fa[b_i + fa.size(0) * 2] * 6000.0 / 6.2831853071795862;
      M[loop_ub] = ma[b_i + ma.size(0) * 2];
      loop_ub = 4 * b_i + 3;
      F[loop_ub] = fa[b_i + fa.size(0) * 3] * 6000.0 / 6.2831853071795862;
      M[loop_ub] = ma[b_i + ma.size(0) * 3];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(i1)

    for (int b_i = 0; b_i < loop_ub; b_i++) {
      F[4 * b_i] = fa[b_i] * 6000.0 / 6.2831853071795862;
      M[4 * b_i] = ma[b_i];
      i1 = 4 * b_i + 1;
      F[i1] = fa[b_i + fa.size(0)] * 6000.0 / 6.2831853071795862;
      M[i1] = ma[b_i + ma.size(0)];
      i1 = 4 * b_i + 2;
      F[i1] = fa[b_i + fa.size(0) * 2] * 6000.0 / 6.2831853071795862;
      M[i1] = ma[b_i + ma.size(0) * 2];
      i1 = 4 * b_i + 3;
      F[i1] = fa[b_i + fa.size(0) * 3] * 6000.0 / 6.2831853071795862;
      M[i1] = ma[b_i + ma.size(0) * 3];
    }
  }
}

//
// File trailer for swsmodel.cpp
//
// [EOF]
//
