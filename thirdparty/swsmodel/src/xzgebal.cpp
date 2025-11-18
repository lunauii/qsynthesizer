//
// File: xzgebal.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "xzgebal.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "coder_array.h"
#include <cmath>
#include <emmintrin.h>

// Function Definitions
//
// Arguments    : double A_data[]
//                int A_size[2]
//                int &ihi
//                double scale_data[]
//                int &scale_size
// Return Type  : int
//
namespace coder {
namespace internal {
namespace reflapack {
int xzgebal(double A_data[], int A_size[2], int &ihi, double scale_data[],
            int &scale_size)
{
  array<double, 2U> b_A_data;
  array<double, 2U> x;
  double temp;
  int b_ix;
  int b_iy;
  int b_k;
  int exitg5;
  int ica;
  int ilo;
  int ix;
  int iy;
  int l;
  int loop_ub;
  int n;
  boolean_T notdone;
  boolean_T skipThisRow;
  n = A_size[0];
  scale_size = n;
  for (int k{0}; k < n; k++) {
    scale_data[k] = 1.0;
  }
  b_k = 0;
  l = n;
  notdone = true;
  do {
    exitg5 = 0;
    if (notdone) {
      int exitg4;
      notdone = false;
      ix = l;
      do {
        exitg4 = 0;
        if (ix > 0) {
          boolean_T exitg6;
          skipThisRow = false;
          ica = 0;
          exitg6 = false;
          while ((!exitg6) && (ica <= static_cast<unsigned char>(l) - 1)) {
            if ((ica + 1 == ix) ||
                (!(A_data[(ix + A_size[0] * ica) - 1] != 0.0))) {
              ica++;
            } else {
              skipThisRow = true;
              exitg6 = true;
            }
          }
          if (skipThisRow) {
            ix--;
          } else {
            scale_data[l - 1] = ix;
            if (ix != l) {
              b_ix = (ix - 1) * n;
              iy = (l - 1) * n;
              ilo = A_size[1];
              x.set_size(n, ilo);
              ica = A_size[0] * A_size[1];
              for (int k{0}; k < ica; k++) {
                x[k] = A_data[k];
              }
              for (int k{0}; k < l; k++) {
                ica = b_ix + k;
                temp = x[ica];
                b_iy = iy + k;
                x[ica] = x[b_iy];
                x[b_iy] = temp;
              }
              for (int k{0}; k < n; k++) {
                ica = k * n;
                b_iy = (ix + ica) - 1;
                temp = x[b_iy];
                ica = (l + ica) - 1;
                x[b_iy] = x[ica];
                x[ica] = temp;
              }
              for (int k{0}; k < ilo; k++) {
                for (int c_k{0}; c_k < n; c_k++) {
                  A_data[c_k + A_size[0] * k] = x[c_k + x.size(0) * k];
                }
              }
            }
            exitg4 = 1;
          }
        } else {
          exitg4 = 2;
        }
      } while (exitg4 == 0);
      if (exitg4 == 1) {
        if (l == 1) {
          ilo = 1;
          ihi = 1;
          exitg5 = 1;
        } else {
          l--;
          notdone = true;
        }
      }
    } else {
      notdone = true;
      while (notdone) {
        boolean_T exitg6;
        notdone = false;
        iy = b_k;
        exitg6 = false;
        while ((!exitg6) && (iy + 1 <= l)) {
          boolean_T exitg7;
          skipThisRow = false;
          ica = b_k;
          exitg7 = false;
          while ((!exitg7) && (ica + 1 <= l)) {
            if ((ica + 1 == iy + 1) ||
                (!(A_data[ica + A_size[0] * iy] != 0.0))) {
              ica++;
            } else {
              skipThisRow = true;
              exitg7 = true;
            }
          }
          if (skipThisRow) {
            iy++;
          } else {
            scale_data[b_k] = iy + 1;
            if (iy + 1 != b_k + 1) {
              b_ix = iy * n;
              ilo = b_k * n;
              loop_ub = A_size[1];
              x.set_size(n, loop_ub);
              ica = A_size[0] * A_size[1];
              for (int c_k{0}; c_k < ica; c_k++) {
                x[c_k] = A_data[c_k];
              }
              for (int c_k{0}; c_k < l; c_k++) {
                ica = b_ix + c_k;
                temp = x[ica];
                b_iy = ilo + c_k;
                x[ica] = x[b_iy];
                x[b_iy] = temp;
              }
              ix = ilo + iy;
              b_iy = ilo + b_k;
              b_ix = n - b_k;
              for (int k{0}; k < b_ix; k++) {
                ica = k * n;
                iy = ix + ica;
                temp = x[iy];
                ica += b_iy;
                x[iy] = x[ica];
                x[ica] = temp;
              }
              for (int k{0}; k < loop_ub; k++) {
                for (int c_k{0}; c_k < n; c_k++) {
                  A_data[c_k + A_size[0] * k] = x[c_k + x.size(0) * k];
                }
              }
            }
            b_k++;
            notdone = true;
            exitg6 = true;
          }
        }
      }
      ilo = b_k + 1;
      ihi = l;
      skipThisRow = false;
      exitg5 = 2;
    }
  } while (exitg5 == 0);
  if (exitg5 != 1) {
    boolean_T exitg3;
    exitg3 = false;
    while ((!exitg3) && (!skipThisRow)) {
      int exitg2;
      skipThisRow = true;
      loop_ub = b_k;
      do {
        exitg2 = 0;
        if (loop_ub + 1 <= l) {
          double c;
          double ca;
          double g;
          double r;
          double s;
          ica = l - b_k;
          ix = loop_ub * n;
          b_A_data.set(&A_data[0], A_size[0], A_size[1]);
          c = blas::xnrm2(ica, b_A_data, (ix + b_k) + 1);
          b_ix = b_k * n + loop_ub;
          iy = b_ix + 1;
          r = 0.0;
          if (ica >= 1) {
            if (ica == 1) {
              r = std::abs(A_data[b_ix]);
            } else {
              temp = 3.3121686421112381E-170;
              ica = (b_ix + (ica - 1) * n) + 1;
              for (int k{iy}; n < 0 ? k >= ica : k <= ica; k += n) {
                g = std::abs(A_data[k - 1]);
                if (g > temp) {
                  s = temp / g;
                  r = r * s * s + 1.0;
                  temp = g;
                } else {
                  s = g / temp;
                  r += s * s;
                }
              }
              r = temp * std::sqrt(r);
              if (std::isnan(r)) {
                b_iy = b_ix + 1;
                int exitg8;
                do {
                  exitg8 = 0;
                  if ((n > 0) && (b_iy <= ica)) {
                    if (std::isnan(A_data[b_iy - 1])) {
                      exitg8 = 1;
                    } else {
                      b_iy += n;
                    }
                  } else {
                    r = rtInf;
                    exitg8 = 1;
                  }
                } while (exitg8 == 0);
              }
            }
          }
          if (l < 1) {
            ica = 0;
          } else {
            ica = 1;
            if (l > 1) {
              temp = std::abs(A_data[ix]);
              for (int k{2}; k <= l; k++) {
                g = std::abs(A_data[(ix + k) - 1]);
                if (g > temp) {
                  ica = k;
                  temp = g;
                }
              }
            }
          }
          ca = std::abs(A_data[(ica + A_size[0] * loop_ub) - 1]);
          ica = n - b_k;
          if (ica < 1) {
            b_iy = 0;
          } else {
            b_iy = 1;
            if (ica > 1) {
              temp = std::abs(A_data[b_ix]);
              for (int k{2}; k <= ica; k++) {
                g = std::abs(A_data[b_ix + (k - 1) * n]);
                if (g > temp) {
                  b_iy = k;
                  temp = g;
                }
              }
            }
          }
          temp = std::abs(A_data[loop_ub + A_size[0] * ((b_iy + b_k) - 1)]);
          if ((c == 0.0) || (r == 0.0)) {
            loop_ub++;
          } else {
            double f;
            int exitg1;
            g = r / 2.0;
            f = 1.0;
            s = c + r;
            do {
              exitg1 = 0;
              if ((c < g) &&
                  (std::fmax(f, std::fmax(c, ca)) < 4.9896007738368E+291) &&
                  (std::fmin(r, std::fmin(g, temp)) >
                   2.0041683600089728E-292)) {
                if (std::isnan(((((c + f) + ca) + r) + g) + temp)) {
                  exitg1 = 1;
                } else {
                  f *= 2.0;
                  c *= 2.0;
                  ca *= 2.0;
                  r /= 2.0;
                  g /= 2.0;
                  temp /= 2.0;
                }
              } else {
                g = c / 2.0;
                while ((g >= r) &&
                       (std::fmax(r, temp) < 4.9896007738368E+291) &&
                       (std::fmin(std::fmin(f, c), std::fmin(g, ca)) >
                        2.0041683600089728E-292)) {
                  f /= 2.0;
                  c /= 2.0;
                  g /= 2.0;
                  ca /= 2.0;
                  r *= 2.0;
                  temp *= 2.0;
                }
                if ((!(c + r >= 0.95 * s)) &&
                    ((!(f < 1.0)) || (!(scale_data[loop_ub] < 1.0)) ||
                     (!(f * scale_data[loop_ub] <= 1.0020841800044864E-292))) &&
                    ((!(f > 1.0)) || (!(scale_data[loop_ub] > 1.0)) ||
                     (!(scale_data[loop_ub] >= 9.9792015476736E+291 / f)))) {
                  temp = 1.0 / f;
                  scale_data[loop_ub] *= f;
                  ica = (b_ix + n * (ica - 1)) + 1;
                  for (int c_k{iy}; n < 0 ? c_k >= ica : c_k <= ica; c_k += n) {
                    A_data[c_k - 1] *= temp;
                  }
                  x.set_size(A_size[0], A_size[1]);
                  ica = A_size[0] * A_size[1];
                  for (int k{0}; k < ica; k++) {
                    x[k] = A_data[k];
                  }
                  b_iy = ix + l;
                  ica = ((((b_iy - ix) / 2) << 1) + ix) + 1;
                  b_ix = ica - 2;
                  for (int k{ix + 1}; k <= b_ix; k += 2) {
                    __m128d b_r;
                    b_r = _mm_loadu_pd(&x[k - 1]);
                    _mm_storeu_pd(&x[k - 1], _mm_mul_pd(_mm_set1_pd(f), b_r));
                  }
                  for (int k{ica}; k <= b_iy; k++) {
                    x[k - 1] = f * x[k - 1];
                  }
                  b_iy = x.size(0);
                  A_size[0] = x.size(0);
                  ica = x.size(1);
                  A_size[1] = x.size(1);
                  for (int k{0}; k < ica; k++) {
                    for (int c_k{0}; c_k < b_iy; c_k++) {
                      b_ix = c_k + x.size(0) * k;
                      A_data[b_ix] = x[b_ix];
                    }
                  }
                  skipThisRow = false;
                }
                exitg1 = 2;
              }
            } while (exitg1 == 0);
            if (exitg1 == 1) {
              exitg2 = 2;
            } else {
              loop_ub++;
            }
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);
      if (exitg2 != 1) {
        exitg3 = true;
      }
    }
  }
  return ilo;
}

} // namespace reflapack
} // namespace internal
} // namespace coder

//
// File trailer for xzgebal.cpp
//
// [EOF]
//
