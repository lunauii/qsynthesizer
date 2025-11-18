//
// File: upfirdnCoreImpl.h
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

#ifndef UPFIRDNCOREIMPL_H
#define UPFIRDNCOREIMPL_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace b_signal {
namespace internal {
namespace upfirdn {
void upfirdnCoreImpl(double xCol, const array<double, 1U> &hCl, double LhD,
                     double pD, double qD, array<double, 1U> &y);

void upfirdnCoreImpl(const array<double, 2U> &xCol,
                     const array<double, 1U> &hCl, double LxD, double LhD,
                     double nChansD, double pD, double qD,
                     array<double, 2U> &y);

} // namespace upfirdn
} // namespace internal
} // namespace b_signal
} // namespace coder

#endif
//
// File trailer for upfirdnCoreImpl.h
//
// [EOF]
//
