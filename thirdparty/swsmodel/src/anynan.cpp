//
// File: anynan.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "anynan.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const double x[8]
// Return Type  : boolean_T
//
namespace coder {
boolean_T anynan(const double x[8])
{
  boolean_T tf;
  tf = false;
  for (int k{0}; k < 8; k++) {
    if (tf || std::isnan(x[k])) {
      tf = true;
    }
  }
  return tf;
}

} // namespace coder

//
// File trailer for anynan.cpp
//
// [EOF]
//
