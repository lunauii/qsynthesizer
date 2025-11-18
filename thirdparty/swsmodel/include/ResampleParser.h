//
// File: ResampleParser.h
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

#ifndef RESAMPLEPARSER_H
#define RESAMPLEPARSER_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
namespace b_signal {
namespace internal {
namespace resample {
class ResampleParser {
public:
  int dim;
  double p;
  double q;
  double inputSize[2];
  boolean_T isRowVectorInput;
  array<double, 1U> filter;
  array<double, 1U> filterWithPadding;
  double filterDelay;
  array<double, 2U> x;
};

} // namespace resample
} // namespace internal
} // namespace b_signal
} // namespace coder

#endif
//
// File trailer for ResampleParser.h
//
// [EOF]
//
