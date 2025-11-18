//
// File: mean.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "mean.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
//
// Arguments    : const array<double, 2U> &x
// Return Type  : double
//
namespace coder {
double mean(const array<double, 2U> &x)
{
  double y;
  if (x.size(1) == 0) {
    y = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x.size(1) <= 1024) {
      firstBlockLength = x.size(1);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(x.size(1)) >> 10);
      lastBlockLength = x.size(1) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x[0];
    for (int k{2}; k <= firstBlockLength; k++) {
      y += x[k - 1];
    }
    for (int k{2}; k <= nblocks; k++) {
      double bsum;
      int hi;
      firstBlockLength = (k - 1) << 10;
      bsum = x[firstBlockLength];
      if (k == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (int b_k{2}; b_k <= hi; b_k++) {
        bsum += x[(firstBlockLength + b_k) - 1];
      }
      y += bsum;
    }
  }
  y /= static_cast<double>(x.size(1));
  return y;
}

} // namespace coder

//
// File trailer for mean.cpp
//
// [EOF]
//
