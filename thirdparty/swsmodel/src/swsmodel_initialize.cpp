//
// File: swsmodel_initialize.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "swsmodel_initialize.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void swsmodel_initialize()
{
  omp_init_nest_lock(&swsmodel_nestLockGlobal);
  isInitialized_swsmodel = true;
}

//
// File trailer for swsmodel_initialize.cpp
//
// [EOF]
//
