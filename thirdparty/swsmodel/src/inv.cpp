//
// File: inv.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "inv.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const double x[64]
//                double y[64]
// Return Type  : void
//
namespace coder {
void inv(const double x[64], double y[64])
{
  double b_x[64];
  double smax;
  int a;
  int ipiv_tmp;
  int jA;
  int kAcol;
  signed char ipiv[8];
  signed char p[8];
  for (int k{0}; k < 64; k++) {
    y[k] = 0.0;
    b_x[k] = x[k];
  }
  for (int k{0}; k < 8; k++) {
    ipiv[k] = static_cast<signed char>(k + 1);
  }
  for (int j{0}; j < 7; j++) {
    int b;
    int jj;
    int jp1j;
    int mmj;
    mmj = 6 - j;
    b = j * 9;
    jj = j * 9;
    jp1j = b + 2;
    jA = 9 - j;
    a = 0;
    smax = std::abs(b_x[jj]);
    for (int k{2}; k < jA; k++) {
      double s;
      s = std::abs(b_x[(b + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (b_x[jj + a] != 0.0) {
      if (a != 0) {
        ipiv_tmp = j + a;
        ipiv[j] = static_cast<signed char>(ipiv_tmp + 1);
        for (int k{0}; k < 8; k++) {
          jA = k << 3;
          kAcol = j + jA;
          smax = b_x[kAcol];
          jA += ipiv_tmp;
          b_x[kAcol] = b_x[jA];
          b_x[jA] = smax;
        }
      }
      jA = (jj - j) + 8;
      for (int k{jp1j}; k <= jA; k++) {
        b_x[k - 1] /= b_x[jj];
      }
    }
    jA = jj;
    for (int b_j{0}; b_j <= mmj; b_j++) {
      smax = b_x[(b + (b_j << 3)) + 8];
      if (smax != 0.0) {
        ipiv_tmp = jA + 10;
        kAcol = (jA - j) + 16;
        for (int ijA{ipiv_tmp}; ijA <= kAcol; ijA++) {
          b_x[ijA - 1] += b_x[((jj + ijA) - jA) - 9] * -smax;
        }
      }
      jA += 8;
    }
  }
  for (int k{0}; k < 8; k++) {
    p[k] = static_cast<signed char>(k + 1);
  }
  for (int k{0}; k < 7; k++) {
    signed char i;
    i = ipiv[k];
    if (i > k + 1) {
      jA = p[i - 1];
      p[i - 1] = p[k];
      p[k] = static_cast<signed char>(jA);
    }
  }
  for (int k{0}; k < 8; k++) {
    jA = (p[k] - 1) << 3;
    y[k + jA] = 1.0;
    for (int b_j{k + 1}; b_j < 9; b_j++) {
      kAcol = (b_j + jA) - 1;
      if (y[kAcol] != 0.0) {
        a = b_j + 1;
        for (int ijA{a}; ijA < 9; ijA++) {
          ipiv_tmp = (ijA + jA) - 1;
          y[ipiv_tmp] -= y[kAcol] * b_x[(ijA + ((b_j - 1) << 3)) - 1];
        }
      }
    }
  }
  for (int k{0}; k < 8; k++) {
    jA = k << 3;
    for (int b_j{7}; b_j >= 0; b_j--) {
      kAcol = b_j << 3;
      a = b_j + jA;
      smax = y[a];
      if (smax != 0.0) {
        y[a] = smax / b_x[b_j + kAcol];
        for (int ijA{0}; ijA < b_j; ijA++) {
          ipiv_tmp = ijA + jA;
          y[ipiv_tmp] -= y[a] * b_x[ijA + kAcol];
        }
      }
    }
  }
}

} // namespace coder

//
// File trailer for inv.cpp
//
// [EOF]
//
