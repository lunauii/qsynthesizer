//
// File: _coder_swsmodel_api.h
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

#ifndef _CODER_SWSMODEL_API_H
#define _CODER_SWSMODEL_API_H

// Include Files
#include "coder_array_mex.h"
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
void swsmodel(coder::array<real_T, 2U> *D, real_T R, real_T H,
              coder::array<real_T, 2U> *F, coder::array<real_T, 2U> *M);

void swsmodel_api(const mxArray *const prhs[3], int32_T nlhs,
                  const mxArray *plhs[2]);

void swsmodel_atexit();

void swsmodel_initialize();

void swsmodel_terminate();

void swsmodel_xil_shutdown();

void swsmodel_xil_terminate();

#endif
//
// File trailer for _coder_swsmodel_api.h
//
// [EOF]
//
