//
// File: xgeev.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "xgeev.h"
#include "rt_nonfinite.h"
#include "xdlahqr.h"
#include "xzgebal.h"
#include "xzlarfg.h"
#include "xzlascl.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const double A_data[]
//                const int A_size[2]
//                creal_T W_data[]
//                int &W_size
// Return Type  : int
//
namespace coder {
namespace internal {
namespace lapack {
int xgeev(const double A_data[], const int A_size[2], creal_T W_data[],
          int &W_size)
{
  array<double, 2U> a;
  array<double, 1U> y;
  double wi_data[8];
  double tau_data[7];
  double absxk;
  double anrm;
  double cfrom1;
  int i;
  int ihi;
  int info;
  int loop_ub;
  int offset;
  boolean_T exitg1;
  loop_ub = A_size[0];
  a.set_size(A_size[0], A_size[1]);
  offset = A_size[0] * A_size[1];
  if (offset - 1 >= 0) {
    std::copy(&A_data[0], &A_data[offset], &a[0]);
  }
  info = 0;
  anrm = 0.0;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i <= offset - 1)) {
    absxk = std::abs(A_data[i]);
    if (std::isnan(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }
      i++;
    }
  }
  if (std::isinf(anrm) || std::isnan(anrm)) {
    W_size = A_size[0];
    for (int j{0}; j < loop_ub; j++) {
      W_data[j].re = rtNaN;
      W_data[j].im = 0.0;
    }
  } else {
    double cscale;
    double ctoc;
    int ilo;
    int n;
    boolean_T guard1;
    boolean_T scalea;
    cscale = anrm;
    scalea = false;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      scalea = true;
      cscale = 6.7178761075670888E-139;
      guard1 = true;
    } else if (anrm > 1.4885657073574029E+138) {
      scalea = true;
      cscale = 1.4885657073574029E+138;
      guard1 = true;
    }
    if (guard1) {
      boolean_T notdone;
      absxk = anrm;
      ctoc = cscale;
      notdone = true;
      while (notdone) {
        double cto1;
        double mul;
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / absxk;
          notdone = false;
        }
        for (int j{0}; j < loop_ub; j++) {
          offset = j * loop_ub - 1;
          for (int b_i{0}; b_i < loop_ub; b_i++) {
            i = (offset + b_i) + 1;
            a[i] = a[i] * mul;
          }
        }
      }
    }
    a.reserve(64);
    y.reserve(8);
    ilo = reflapack::xzgebal((double *)a.data(), a.size(), ihi,
                             (double *)y.data(), offset);
    y.set_size(offset);
    n = a.size(0);
    if ((ihi - ilo) + 1 > 1) {
      offset = static_cast<unsigned char>(ilo - 1);
      if (offset - 1 >= 0) {
        std::memset(&tau_data[0], 0,
                    static_cast<unsigned int>(offset) * sizeof(double));
      }
      if (ihi <= n - 1) {
        std::memset(&tau_data[ihi + -1], 0,
                    static_cast<unsigned int>(n - ihi) * sizeof(double));
      }
      offset = static_cast<signed char>(a.size(0));
      y.set_size(static_cast<int>(static_cast<signed char>(a.size(0))));
      if (offset - 1 >= 0) {
        std::memset(&y[0], 0,
                    static_cast<unsigned int>(offset) * sizeof(double));
      }
      for (int c_i{ilo}; c_i < ihi; c_i++) {
        int b_lastv;
        int d_i;
        int exitg2;
        int in;
        int iv0;
        int jA;
        int lastc;
        int lastv;
        offset = (c_i - 1) * n;
        in = c_i * n;
        lastv = ihi - c_i;
        cfrom1 = a[c_i + a.size(0) * (c_i - 1)];
        i = c_i + 2;
        if (i > n) {
          i = n;
        }
        absxk = reflapack::xzlarfg(lastv, cfrom1, a, i + offset);
        tau_data[c_i - 1] = absxk;
        d_i = c_i + a.size(0) * (c_i - 1);
        a[d_i] = 1.0;
        iv0 = c_i + offset;
        jA = in + 1;
        if (absxk != 0.0) {
          b_lastv = lastv;
          offset = iv0 + lastv;
          while ((b_lastv > 0) && (a[offset - 1] == 0.0)) {
            b_lastv--;
            offset--;
          }
          lastc = ihi;
          exitg1 = false;
          while ((!exitg1) && (lastc > 0)) {
            offset = in + lastc;
            i = offset;
            do {
              exitg2 = 0;
              if ((n > 0) && (i <= offset + (b_lastv - 1) * n)) {
                if (a[i - 1] != 0.0) {
                  exitg2 = 1;
                } else {
                  i += n;
                }
              } else {
                lastc--;
                exitg2 = 2;
              }
            } while (exitg2 == 0);
            if (exitg2 == 1) {
              exitg1 = true;
            }
          }
        } else {
          b_lastv = 0;
          lastc = 0;
        }
        if (b_lastv > 0) {
          if (lastc != 0) {
            offset = static_cast<unsigned char>(lastc);
            std::memset(&y[0], 0,
                        static_cast<unsigned int>(offset) * sizeof(double));
            offset = iv0;
            i = (in + n * (b_lastv - 1)) + 1;
            for (int j{jA}; n < 0 ? j >= i : j <= i; j += n) {
              loop_ub = j + lastc;
              for (int b_i{j}; b_i < loop_ub; b_i++) {
                info = b_i - j;
                y[info] = y[info] + a[b_i - 1] * a[offset];
              }
              offset++;
            }
          }
          ctoc = -tau_data[c_i - 1];
          if (!(ctoc == 0.0)) {
            info = in;
            offset = static_cast<unsigned char>(b_lastv);
            for (int j{0}; j < offset; j++) {
              absxk = a[iv0 + j];
              if (absxk != 0.0) {
                absxk *= ctoc;
                i = info + 1;
                loop_ub = lastc + info;
                for (int b_i{i}; b_i <= loop_ub; b_i++) {
                  a[b_i - 1] = a[b_i - 1] + y[(b_i - info) - 1] * absxk;
                }
              }
              info += n;
            }
          }
        }
        jA = (c_i + in) + 1;
        ctoc = tau_data[c_i - 1];
        if (ctoc != 0.0) {
          i = iv0 + lastv;
          while ((lastv > 0) && (a[i - 1] == 0.0)) {
            lastv--;
            i--;
          }
          info = (n - c_i) - 1;
          exitg1 = false;
          while ((!exitg1) && (info + 1 > 0)) {
            offset = jA + info * n;
            i = offset;
            do {
              exitg2 = 0;
              if (i <= (offset + lastv) - 1) {
                if (a[i - 1] != 0.0) {
                  exitg2 = 1;
                } else {
                  i++;
                }
              } else {
                info--;
                exitg2 = 2;
              }
            } while (exitg2 == 0);
            if (exitg2 == 1) {
              exitg1 = true;
            }
          }
        } else {
          lastv = 0;
          info = -1;
        }
        if (lastv > 0) {
          if (info + 1 != 0) {
            if (info >= 0) {
              std::memset(&y[0], 0,
                          static_cast<unsigned int>(info + 1) * sizeof(double));
            }
            offset = 0;
            i = jA + n * info;
            for (int j{jA}; n < 0 ? j >= i : j <= i; j += n) {
              absxk = 0.0;
              loop_ub = j + lastv;
              for (int b_i{j}; b_i < loop_ub; b_i++) {
                absxk += a[b_i - 1] * a[(iv0 + b_i) - j];
              }
              y[offset] = y[offset] + absxk;
              offset++;
            }
          }
          if (!(-ctoc == 0.0)) {
            for (int j{0}; j <= info; j++) {
              if (y[j] != 0.0) {
                absxk = y[j] * -ctoc;
                offset = lastv + jA;
                for (int b_i{jA}; b_i < offset; b_i++) {
                  a[b_i - 1] = a[b_i - 1] + a[(iv0 + b_i) - jA] * absxk;
                }
              }
              jA += n;
            }
          }
        }
        a[d_i] = cfrom1;
      }
    }
    a.reserve(64);
    y.reserve(8);
    info = reflapack::xdlahqr(ilo, ihi, (double *)a.data(), a.size(),
                              (double *)y.data(), i, wi_data, offset);
    y.set_size(i);
    if (scalea) {
      offset = A_size[0] - info;
      y.reserve(8);
      reflapack::xzlascl(cscale, anrm, offset, (double *)y.data(), info + 1);
      reflapack::xzlascl(cscale, anrm, offset, wi_data, info + 1);
      if (info != 0) {
        y.reserve(8);
        reflapack::xzlascl(cscale, anrm, ilo - 1, (double *)y.data(), 1);
        reflapack::xzlascl(cscale, anrm, ilo - 1, wi_data, 1);
      }
    }
    if (info != 0) {
      for (int j{ilo}; j <= info; j++) {
        y[j - 1] = rtNaN;
        wi_data[j - 1] = 0.0;
      }
    }
    offset = y.size(0);
    W_size = y.size(0);
    for (int j{0}; j < offset; j++) {
      W_data[j].re = y[j];
      W_data[j].im = wi_data[j];
    }
  }
  return info;
}

} // namespace lapack
} // namespace internal
} // namespace coder

//
// File trailer for xgeev.cpp
//
// [EOF]
//
