#include "poly_coeff.h"

/*****************************************************************
MODULE:   get_3rd_order_poly_coeff
*****************************************************************/
void get_3rd_order_poly_coeff
(
    float aot[], /* I: array of AOT values [NAOT_VALS] */
    float *atm,  /* I: atmospheric variable, corresponding to aerosol vector */
    int nAtm,    /* I: dimension of atm array, max NAOT_VALS  */
    float *coeff /* O: regression coefficient, quadratic [NCOEF] */
)
{
    float y1[NCOEF] = {0.0, 0.0, 0.0, 0.0};
    float x[NAOT_VALS][NCOEF];
    float z[NCOEF][NCOEF], z1[NCOEF*NCOEF], zinv[NCOEF][NCOEF], zinv1[16];
    int i,j,k;
  
    for (i=0; i<nAtm; i++)
    {
        x[i][0] = aot[i] * aot[i] * aot[i];
        x[i][1] = aot[i] * aot[i];
        x[i][2] = aot[i];
        x[i][3] = 1.0;
    }
  
    for (i=0; i<NCOEF; i++)
        for (j=0; j<NCOEF; j++)
        {
           z[i][j] = 0.0;
           for (k=0; k<nAtm; k++)
               z[i][j] += x[k][i] * x[k][j];
           z1[NCOEF*i + j] = z[i][j];
        }
  
    inverseMatrix4x4 (z1, zinv1);
  
    for (i=0; i<NCOEF; i++)
        for (j=0; j<NCOEF; j++)
            zinv[i][j] = zinv1[NCOEF*i + j];
  
    for (i=0; i<NCOEF; i++)
        for (k=0; k<nAtm; k++)
            y1[i] += x[k][i] * atm[k];
  
    for (i=0; i<NCOEF; i++)
    {
        coeff[i] = 0.0;
        for (j=0; j<NCOEF; j++)
            coeff[i] += zinv[i][j] * y1[j];
    }
}


/*****************************************************************
MODULE:   invf
*****************************************************************/
float invf
(
    int i,
    int j,
    const float *m
)
{
    int o = 2 + (j-i);
    i += 4+o;
    j += 4-o;
  
    #define e(a,b) m[ ((j+b)%4)*4 + ((i+a)%4) ]
  
    float inv =
      e(+1,-1)*e(+0,+0)*e(-1,+1) +
      e(+1,+1)*e(+0,-1)*e(-1,+0) +
      e(-1,-1)*e(+1,+0)*e(+0,+1) -
      e(-1,-1)*e(+0,+0)*e(+1,+1) -
      e(-1,+1)*e(+0,-1)*e(+1,+0) -
      e(+1,-1)*e(-1,+0)*e(+0,+1);
  
    return (o%2)?inv : -inv;
  
    #undef e
}


/*****************************************************************
MODULE:   inverseMatrix4x4
*****************************************************************/
bool inverseMatrix4x4
(
    const float *m,
    float *out
)
{
    int i,j,k;
    float inv[NCOEF*NCOEF];
    float D = 0;

    for (i=0; i<NCOEF; i++)
        for (j=0; j<NCOEF; j++)
            inv[j*NCOEF+i] = invf(i,j,m);
  
    for (k=0; k<NCOEF; k++)
        D += m[k] * inv[k*NCOEF];
  
    if (D == 0)
        return false;
  
    D = 1.0 / D;
  
    for (i=0; i<NCOEF*NCOEF; i++)
        out[i] = inv[i] * D;
  
    return true;
}

