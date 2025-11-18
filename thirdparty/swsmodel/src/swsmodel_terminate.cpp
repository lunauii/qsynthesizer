//
// File: swsmodel_terminate.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "swsmodel_terminate.h"
#include "rt_nonfinite.h"
#include "swsmodel_data.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void swsmodel_terminate()
{
  omp_destroy_nest_lock(&swsmodel_nestLockGlobal);
  isInitialized_swsmodel = false;
}

//
// File trailer for swsmodel_terminate.cpp
//
// [EOF]
//
