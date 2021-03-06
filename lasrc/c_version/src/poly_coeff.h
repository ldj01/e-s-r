#ifndef _POLY_COEFF_H_
#define _POLY_COEFF_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "common.h"

/* Prototypes */
void get_3rd_order_poly_coeff (float aot[], float *atm, int nAtm, float *coeff);

float invf (int i, int j, const float* m);

bool inverseMatrix4x4 (const float *m, float *out);

#endif
