//
// File: filter.h
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

#ifndef FILTER_H
#define FILTER_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void filter(const array<double, 2U> &x, array<double, 2U> &y);

void filter(const double b[9], const array<double, 2U> &x,
            array<double, 2U> &y);

} // namespace coder

#endif
//
// File trailer for filter.h
//
// [EOF]
//
