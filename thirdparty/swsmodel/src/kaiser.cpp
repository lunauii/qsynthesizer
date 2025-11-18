//
// File: kaiser.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "kaiser.h"
#include "casyi.h"
#include "cmlri.h"
#include "gammaln.h"
#include "log.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// Arguments    : double N
//                array<double, 1U> &w
// Return Type  : void
//
namespace coder {
void kaiser(double N, array<double, 1U> &w)
{
  creal_T tmp;
  creal_T zd;
  double r;
  int inw;
  int nw;
  if (N == std::floor(N)) {
    nw = static_cast<int>(N);
  } else {
    nw = static_cast<int>(std::round(N));
  }
  w.set_size(nw);
  if (nw <= 1) {
    w.set_size(nw);
    for (int k{0}; k < nw; k++) {
      w[k] = 1.0;
    }
  } else {
    int iseven;
    int mid;
    iseven = 1 - static_cast<int>(static_cast<unsigned int>(nw) & 1U);
    mid = (nw >> 1) + 1;
    if (mid <= nw) {
      zd.im = 0.0;
    }
    for (int k{mid}; k <= nw; k++) {
      double ak;
      double s;
      double zd_tmp;
      r = static_cast<double>(iseven + ((k - mid) << 1)) /
          (static_cast<double>(nw) - 1.0);
      zd_tmp = 5.0 * std::sqrt((1.0 - r) * (r + 1.0));
      zd.re = zd_tmp;
      if (std::isnan(zd_tmp)) {
        tmp.re = rtNaN;
        tmp.im = 0.0;
      } else {
        double az;
        int b_nw;
        int ierr;
        boolean_T guard1;
        ierr = 0;
        if (zd_tmp > 0.0) {
          r = zd_tmp;
        } else {
          r = 0.0;
        }
        if (r > 1.0737418235E+9) {
          ierr = 4;
        } else if (r > 32767.999992370605) {
          ierr = 3;
        }
        tmp.re = 0.0;
        tmp.im = 0.0;
        if (zd_tmp > 0.0) {
          az = zd_tmp;
        } else {
          az = 0.0;
        }
        guard1 = false;
        if (az <= 2.0) {
          b_nw = 0;
          if (zd_tmp > 0.0) {
            if (zd_tmp < 2.2250738585072014E-305) {
              tmp.re = 1.0;
              tmp.im = 0.0;
            } else {
              double acz;
              double cz_re;
              tmp.re = 0.5 * zd_tmp;
              tmp.im = 0.0;
              if (zd_tmp > 4.7170688552396617E-153) {
                cz_re = tmp.re * tmp.re;
                if (cz_re > 0.0) {
                  acz = cz_re;
                } else {
                  acz = 0.0;
                }
              } else {
                cz_re = 0.0;
                acz = 0.0;
              }
              b_log(tmp);
              r = 1.0;
              gammaln(r);
              r = tmp.re * 0.0 - r;
              ak = tmp.im * 0.0;
              if (r > -700.92179369444591) {
                double coef_im;
                double coef_re;
                double s1_im;
                double s1_re;
                r = std::exp(r);
                coef_re = r * std::cos(ak);
                coef_im = r * std::sin(ak);
                r = 2.2204460492503131E-16 * acz;
                s1_re = 1.0;
                s1_im = 0.0;
                if (!(acz < 2.2204460492503131E-16)) {
                  double aa;
                  tmp.re = 1.0;
                  tmp.im = 0.0;
                  ak = 3.0;
                  s = 1.0;
                  aa = 2.0;
                  double im;
                  double re;
                  double rs;
                  do {
                    rs = 1.0 / s;
                    re = tmp.re * cz_re - tmp.im * 0.0;
                    im = tmp.re * 0.0 + tmp.im * cz_re;
                    tmp.re = rs * re;
                    tmp.im = rs * im;
                    s1_re += tmp.re;
                    s1_im += tmp.im;
                    s += ak;
                    ak += 2.0;
                    aa = aa * acz * rs;
                  } while (!!(aa > r));
                }
                ak = s1_re * coef_re - s1_im * coef_im;
                r = s1_re * coef_im + s1_im * coef_re;
                tmp.re = ak - r * 0.0;
                tmp.im = ak * 0.0 + r;
              } else {
                b_nw = 1;
                tmp.re = 0.0;
                tmp.im = 0.0;
                if (acz > 0.0) {
                  b_nw = -1;
                }
              }
            }
          } else {
            tmp.re = 1.0;
            tmp.im = 0.0;
          }
          if (b_nw < 0) {
            inw = 1;
          } else {
            inw = b_nw;
          }
          if ((1 - inw != 0) && (b_nw < 0)) {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
        if (guard1) {
          if (az < 21.784271729432426) {
            b_nw = cmlri(zd, tmp);
            if (b_nw < 0) {
              if (b_nw == -2) {
                inw = -2;
              } else {
                inw = -1;
              }
            } else {
              inw = 0;
            }
          } else {
            b_nw = casyi(zd, tmp);
            if (b_nw < 0) {
              if (b_nw == -2) {
                inw = -2;
              } else {
                inw = -1;
              }
            } else {
              inw = 0;
            }
          }
        }
        guard1 = false;
        if (inw < 0) {
          if (inw == -2) {
            tmp.re = rtNaN;
            tmp.im = 0.0;
          } else {
            ierr = 2;
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
        if (guard1 && (ierr == 2)) {
          tmp.re = rtInf;
          tmp.im = 0.0;
        }
        if (zd_tmp > 0.0) {
          r = tmp.re;
          tmp.re = r;
          tmp.im = 0.0;
        }
      }
      if (tmp.im == 0.0) {
        r = tmp.re / 27.239871823604449;
        ak = 0.0;
      } else if (tmp.re == 0.0) {
        r = 0.0;
        ak = tmp.im / 27.239871823604449;
      } else {
        r = tmp.re / 27.239871823604449;
        ak = tmp.im / 27.239871823604449;
      }
      s = std::abs(r);
      r = std::abs(ak);
      if (s < r) {
        s /= r;
        w[k - 1] = r * std::sqrt(s * s + 1.0);
      } else if (s > r) {
        r /= s;
        w[k - 1] = s * std::sqrt(r * r + 1.0);
      } else if (std::isnan(r)) {
        w[k - 1] = rtNaN;
      } else {
        w[k - 1] = s * 1.4142135623730951;
      }
    }
    for (int k{0}; k <= mid - 2; k++) {
      w[k] = w[(nw - k) - 1];
    }
  }
}

} // namespace coder

//
// File trailer for kaiser.cpp
//
// [EOF]
//
