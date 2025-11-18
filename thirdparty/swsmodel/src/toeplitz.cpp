//
// File: toeplitz.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "toeplitz.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : const double c[8]
//                double t[64]
// Return Type  : void
//
namespace coder {
void toeplitz(const double c[8], double t[64])
{
  int ij;
  ij = 0;
  for (int j{0}; j < 8; j++) {
    int k;
    k = j;
    for (int i{0}; i < 8; i++) {
      if (i < j) {
        t[ij + i] = c[k];
        k--;
      } else {
        t[ij + i] = c[k];
        k++;
      }
    }
    ij += 8;
  }
}

} // namespace coder

//
// File trailer for toeplitz.cpp
//
// [EOF]
//
