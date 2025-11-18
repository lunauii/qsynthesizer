//
// File: uniformResampleKernel.h
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

#ifndef UNIFORMRESAMPLEKERNEL_H
#define UNIFORMRESAMPLEKERNEL_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace coder {
namespace b_signal {
namespace internal {
namespace resample {
class ResampleParser;

}
} // namespace internal
} // namespace b_signal
} // namespace coder

// Function Declarations
namespace coder {
namespace b_signal {
namespace internal {
namespace resample {
void uniformResampleAlongFirstDim(array<double, 2U> &xIn,
                                  const ResampleParser &opts);

}
} // namespace internal
} // namespace b_signal
} // namespace coder

#endif
//
// File trailer for uniformResampleKernel.h
//
// [EOF]
//
