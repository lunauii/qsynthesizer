//
// File: mldivide.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "mldivide.h"
#include "qrsolve.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const array<double, 2U> &A
//                array<double, 1U> &B
// Return Type  : void
//
namespace coder {
void mldivide(const array<double, 2U> &A, array<double, 1U> &B)
{
  array<double, 2U> b_A;
  array<double, 1U> b_B;
  array<int, 2U> ipiv;
  if ((A.size(0) == 0) || (A.size(1) == 0) || (B.size(0) == 0)) {
    int yk;
    yk = A.size(1);
    B.set_size(A.size(1));
    for (int k{0}; k < yk; k++) {
      B[k] = 0.0;
    }
  } else if (A.size(0) == A.size(1)) {
    double temp;
    int LDA;
    int jA;
    int n;
    int yk;
    yk = A.size(0);
    n = A.size(1);
    if (yk <= n) {
      n = yk;
    }
    yk = B.size(0);
    if (yk <= n) {
      n = yk;
    }
    LDA = A.size(0);
    b_A.set_size(A.size(0), A.size(1));
    yk = A.size(0) * A.size(1);
    for (int k{0}; k < yk; k++) {
      b_A[k] = A[k];
    }
    ipiv.set_size(1, n);
    if (n > 0) {
      ipiv[0] = 1;
      yk = 1;
      for (int k{2}; k <= n; k++) {
        yk++;
        ipiv[k - 1] = yk;
      }
    }
    if (n >= 1) {
      int u0;
      u0 = n - 1;
      if (u0 > n) {
        u0 = n;
      }
      for (int j{0}; j < u0; j++) {
        int b;
        int jp1j;
        int mmj;
        int temp_tmp;
        mmj = n - j;
        b = j * (LDA + 1);
        jp1j = b + 2;
        if (mmj - 1 < 0) {
          yk = -1;
        } else {
          yk = 0;
          if (mmj - 1 > 0) {
            temp = std::abs(b_A[b]);
            for (int k{2}; k <= mmj; k++) {
              double s;
              s = std::abs(b_A[(b + k) - 1]);
              if (s > temp) {
                yk = k - 1;
                temp = s;
              }
            }
          }
        }
        if (b_A[b + yk] != 0.0) {
          if (yk != 0) {
            jA = j + yk;
            ipiv[j] = jA + 1;
            for (int k{0}; k < n; k++) {
              yk = k * LDA;
              temp_tmp = j + yk;
              temp = b_A[temp_tmp];
              yk += jA;
              b_A[temp_tmp] = b_A[yk];
              b_A[yk] = temp;
            }
          }
          yk = b + mmj;
          for (int k{jp1j}; k <= yk; k++) {
            b_A[k - 1] = b_A[k - 1] / b_A[b];
          }
        }
        yk = b + LDA;
        jA = yk;
        for (int k{0}; k <= mmj - 2; k++) {
          temp = b_A[yk + k * LDA];
          if (temp != 0.0) {
            temp_tmp = jA + 2;
            jp1j = mmj + jA;
            for (int ijA{temp_tmp}; ijA <= jp1j; ijA++) {
              b_A[ijA - 1] = b_A[ijA - 1] + b_A[((b + ijA) - jA) - 1] * -temp;
            }
          }
          jA += LDA;
        }
      }
    }
    for (int k{0}; k <= n - 2; k++) {
      yk = ipiv[k];
      if (yk != k + 1) {
        temp = B[k];
        B[k] = B[yk - 1];
        B[yk - 1] = temp;
      }
    }
    for (int k{0}; k < n; k++) {
      yk = LDA * k;
      if (B[k] != 0.0) {
        jA = k + 2;
        for (int ijA{jA}; ijA <= n; ijA++) {
          B[ijA - 1] = B[ijA - 1] - B[k] * b_A[(ijA + yk) - 1];
        }
      }
    }
    for (int k{n}; k >= 1; k--) {
      yk = LDA * (k - 1);
      temp = B[k - 1];
      if (temp != 0.0) {
        temp /= b_A[(k + yk) - 1];
        B[k - 1] = temp;
        for (int ijA{0}; ijA <= k - 2; ijA++) {
          B[ijA] = B[ijA] - B[k - 1] * b_A[ijA + yk];
        }
      }
    }
  } else {
    int yk;
    b_B.set_size(B.size(0));
    yk = B.size(0) - 1;
    for (int k{0}; k <= yk; k++) {
      b_B[k] = B[k];
    }
    internal::qrsolve(A, b_B, B);
  }
}

} // namespace coder

//
// File trailer for mldivide.cpp
//
// [EOF]
//
