//
// File: roots.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "roots.h"
#include "rt_nonfinite.h"
#include "xgeev.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const double c[9]
//                creal_T r_data[]
// Return Type  : int
//
namespace coder {
int roots(const double c[9], creal_T r_data[])
{
  creal_T eiga_data[8];
  double a_data[64];
  int a_size[2];
  int k2;
  int nTrailingZeros;
  int r_size;
  std::memset(&r_data[0], 0, 8U * sizeof(creal_T));
  r_size = 1;
  while ((r_size <= 9) && (!(c[r_size - 1] != 0.0))) {
    r_size++;
  }
  k2 = 9;
  while ((k2 >= r_size) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }
  nTrailingZeros = 9 - k2;
  if (r_size < k2) {
    double ctmp[9];
    int companDim;
    boolean_T exitg1;
    companDim = k2 - r_size;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      int j;
      boolean_T exitg2;
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[r_size + j] / c[r_size - 1];
        if (std::isinf(std::abs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }
      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        r_size++;
        companDim--;
      }
    }
    if (companDim < 1) {
      r_size = 9 - k2;
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      r_size = companDim * companDim;
      std::memset(&a_data[0], 0,
                  static_cast<unsigned int>(r_size) * sizeof(double));
      for (int k{0}; k <= companDim - 2; k++) {
        r_size = companDim * k;
        a_data[r_size] = -ctmp[k];
        a_data[(k + r_size) + 1] = 1.0;
      }
      a_data[companDim * (companDim - 1)] = -ctmp[companDim - 1];
      if (nTrailingZeros - 1 >= 0) {
        std::memset(&r_data[0], 0,
                    static_cast<unsigned int>(nTrailingZeros) *
                        sizeof(creal_T));
      }
      if (companDim == 1) {
        for (int k{0}; k < companDim; k++) {
          eiga_data[k].re = a_data[k];
          eiga_data[k].im = 0.0;
        }
      } else {
        internal::lapack::xgeev(a_data, a_size, eiga_data, r_size);
      }
      for (int k{0}; k < companDim; k++) {
        r_data[(k - k2) + 9] = eiga_data[k];
      }
      r_size = (companDim - k2) + 9;
    }
  } else {
    r_size = 9 - k2;
  }
  return r_size;
}

} // namespace coder

//
// File trailer for roots.cpp
//
// [EOF]
//
