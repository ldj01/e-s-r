#include "ar.h"
#include "const.h"
#include "error.h"
#include "sixs_runs.h"
#include "clouds.h"

/* #define VRA_THRESHOLD 0.1 */
#define VRA_THRESHOLD 0.08

extern atmos_t atmos_coef;

int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef);
int free_mem_atmos_coeff(atmos_t *atmos_coef);
void SrInterpAtmCoef (Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef, atmos_t *interpol_atmos_coef);


bool cloud_detection_pass1
(
    Lut_t *lut,              /* I: lookup table informat */
    int nsamp,               /* I: number of samples to be processed */
    int il,                  /* I: current line being processed */
    uint16_t **line_in,      /* I: array of input lines, one for each band */
    uint8 *qa_line,          /* I: array of QA data for the current line */
    uint16_t *b6_line,       /* I: array of thermal data for the current line */
    float *atemp_line,       /* I: auxiliary temperature for the line */
    cld_diags_t *cld_diags   /* I/O: cloud diagnostics (stats are updated) */
)
{
    int is;                   /* current sample in the line */
    bool is_fill;             /* is the current pixel fill */
    float tmpflt;             /* temporary floating point value */
    float rho1, rho3, rho4, rho5, rho7, t6;  /* reflectance and temp values */
    int C1, C2=0, C3, C4=0, C5=0, water;  /* cloud and water
                                                            indicators */
    int cld_row, cld_col;     /* cloud line, sample location */
    float vra, ndvi;          /* NDVI value */
    Img_coord_int_t loc;      /* line/sample location for current pix */
    atmos_t interpol_atmos_coef; /* interpolated atmospheric coefficients,
                                    based on the current line/sample location
                                    in the aerosol data grid */

    /* Allocate memory for the interpolated atmospheric coefficients and
       start the location for the current line and current cloud row */
    allocate_mem_atmos_coeff (1, &interpol_atmos_coef);
    loc.l = il;
    cld_row = il / cld_diags->cellheight;

    int thresh_TM;
    /* Calculate a threshold that was a hardcoded magic number before.
       This was set to 5000 with the original scaling factor. */
    thresh_TM = 0.5 * lut->mult_factor + lut->add_offset;

    /* Loop through the samples in this line */
    for (is = 0; is < nsamp; is++) {
        loc.s = is;
        cld_col = is / cld_diags->cellwidth;

        if ((qa_line[is]&0x01)==0x01)
            is_fill=true;
        else
            is_fill=false;

        if (!is_fill) {
            if (((qa_line[is] & 0x08) == 0x00) ||
              ((lut->meta.inst == INST_TM) && (line_in[2][is] < thresh_TM)))
            { /* no saturation in band 3 */
                /* Interpolate the atmospheric coefficients for the current
                   pixel */
                SrInterpAtmCoef(lut, &loc, &atmos_coef, &interpol_atmos_coef);

                /* Compute the reflectance for each band using the interpolated
                   atmospheric coefficients */
                rho1 = line_in[0][is] * lut->scale_factor + lut->add_offset;
                rho1 = (rho1/interpol_atmos_coef.tgOG[0][0] -
                    interpol_atmos_coef.rho_ra[0][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[0][0] *
                    interpol_atmos_coef.td_ra[0][0] *
                    interpol_atmos_coef.tu_ra[0][0]);
                rho1 /= tmpflt;
                rho1 /= (1. + interpol_atmos_coef.S_ra[0][0] * rho1);

                rho3 = line_in[2][is] * lut->scale_factor + lut->add_offset;
                rho3 = (rho3 / interpol_atmos_coef.tgOG[2][0] -
                    interpol_atmos_coef.rho_ra[2][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[2][0] *
                    interpol_atmos_coef.td_ra[2][0] *
                    interpol_atmos_coef.tu_ra[2][0]);
                rho3 /= tmpflt;
                rho3 /= (1. + interpol_atmos_coef.S_ra[2][0] * rho3);

                rho4 = line_in[3][is] * lut->scale_factor + lut->add_offset;
                rho4 = (rho4 / interpol_atmos_coef.tgOG[3][0] -
                    interpol_atmos_coef.rho_ra[3][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[3][0] *
                    interpol_atmos_coef.td_ra[3][0] *
                    interpol_atmos_coef.tu_ra[3][0]);
                rho4 /= tmpflt;
                rho4 /= (1. + interpol_atmos_coef.S_ra[3][0] * rho4);

                rho5 = line_in[4][is] * lut->scale_factor + lut->add_offset;
                rho5 = (rho5 / interpol_atmos_coef.tgOG[4][0] -
                    interpol_atmos_coef.rho_ra[4][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[4][0] *
                    interpol_atmos_coef.td_ra[4][0] *
                    interpol_atmos_coef.tu_ra[4][0]);
                rho5 /= tmpflt;
                rho5 /= (1. + interpol_atmos_coef.S_ra[4][0] * rho5);

                rho7 = line_in[5][is] * lut->scale_factor + lut->add_offset;
                rho7 = (rho7 / interpol_atmos_coef.tgOG[5][0] -
                    interpol_atmos_coef.rho_ra[5][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[5][0] *
                    interpol_atmos_coef.td_ra[5][0] *
                    interpol_atmos_coef.tu_ra[5][0]);
                rho7 /= tmpflt;
                rho7 /= (1. + interpol_atmos_coef.S_ra[5][0] * rho7);

                /* Get the temperature */
                t6 = b6_line[is] * lut->b6_scale_factor + lut->b6_add_offset;

                /* Compute cloud coefficients */
                vra = rho1 - rho3 * 0.5;
                    
                C1 = (int)(vra > VRA_THRESHOLD);
                C2 = (t6 < (atemp_line[is]-7.));  
   
                tmpflt = rho4 / rho3;
                C3 = ((tmpflt >= 0.9) && (tmpflt <= 1.3));

                C4 = (rho7 > 0.03);
                C5 = ((rho3 > 0.6)||(rho4 > 0.6));
                    
                /**
                Water test :
                ndvi < 0 => water
                ((0<ndvi<0.1) or (b4<5%)) and b5 < 0.01 => turbid water
                **/
                if (rho4 + rho3 != 0)
                    ndvi = (rho4 - rho3) / (rho4 + rho3);
                else
                    ndvi = 0.01;
                water = (ndvi < 0) || ((((ndvi > 0) && (ndvi < 0.1)) ||
                    (rho4 < 0.05)) && (rho5 < 0.02));
                
                if (!water) { /* if not water */
                    if ((t6 > (atemp_line[is] - 20.)) && (!C5)) { 
                        if (!((C1||C3)&&C2&&C4)) { /* clear */
                            cld_diags->avg_t6_clear[cld_row][cld_col] += t6;
                            cld_diags->std_t6_clear[cld_row][cld_col] +=
                                (t6*t6);
                            cld_diags->avg_b7_clear[cld_row][cld_col] += rho7;
                            cld_diags->std_b7_clear[cld_row][cld_col] +=
                                (rho7*rho7);
                            cld_diags->nb_t6_clear[cld_row][cld_col]++;
                        }
                    }
                }
            }  /* end if no saturation in band 3 */
        }  /* end if !is_fill */ 
    }  /* end for is */

    free_mem_atmos_coeff(&interpol_atmos_coef);
    return true;
}


bool cloud_detection_pass2
(
    Lut_t *lut,              /* I: lookup table informat */
    int nsamp,               /* I: number of samples to be processed */
    int il,                  /* I: current line being processed */
    uint16_t **line_in,      /* I: array of input lines, one for each band */
    uint8 *qa_line,          /* I: array of QA data for the current line */
    uint16_t *b6_line,       /* I: array of thermal data for the current line */
    cld_diags_t *cld_diags,  /* I: cloud diagnostics */
    char *ddv_line           /* O: dark dense vegetation line */
            /**
            use ddv_line to store internal cloud screening info
            bit 2 = adjacent cloud 1=yes 0=no
            bit 3 = fill value 1=fill 0=valid
            bit 4 = land/water mask 1=land 0=water
            bit 5 = cloud 0=clear 1=cloudy
            bit 6 = cloud shadow 
            bit 7 = snow
            **/
)
{
    int is;                   /* current sample in the line */
    bool is_fill;             /* is the current pixel fill */
    int il_ar, is_ar;         /* line/sample in the aerosol region */
    bool thermal_band;        /* is thermal data available */
    float rho1, rho2, rho3, rho4, rho5, rho7;  /* reflectance & temp values */
    float t6 = 0.0;           /* temperature values */
    int C1, C2=0, C4=0, C5=0, water;  /* cloud and water
                                                            indicators */
    int cld_row, cld_col;     /* cloud line, sample location */
    float vra, ndvi, ndsi, temp_snow_thshld; /* NDVI, NDSI, snow threshold */
    float temp_b6_clear,temp_thshld1,temp_thshld2,atemp_ancillary;
    float tmpflt, tmpflt_arr[10];  /* temporary floats */
    Img_coord_int_t loc;      /* line/sample location for current pix */
    atmos_t interpol_atmos_coef; /* interpolated atmospheric coefficients,
                                    based on the current line/sample location
                                    in the aerosol data grid */
    int band5_thresh;
    int thresh_TM;

    /* Initialize the thermal band information and snow threshold */
    thermal_band = true;
    if (b6_line == NULL) 
        thermal_band = false;

    /* This appears to be an unscaled value which is outside the
       maximum value for this band - DLE */
    temp_snow_thshld = 380.;  /* now flag snow and possibly salt pan */

    /* Calculate a threshold that was a hardcoded magic number before.
       This was set to 2000 with the original scaling factor. */
    band5_thresh = 0.2 * lut->mult_factor + lut->add_offset;

    /* Calculate a threshold that was a hardcoded magic number before.
       This was set to 5000 with the original scaling factor. */
    thresh_TM = 0.5 * lut->mult_factor + lut->add_offset;

    /* Allocate memory for the interpolated atmospheric coefficients and
       start the location for the current line and current cloud row */
    allocate_mem_atmos_coeff(1,&interpol_atmos_coef);
    loc.l = il;
    cld_row = il / cld_diags->cellheight;

    /* Initialize the aerosol line location value */
    il_ar = il / lut->ar_region_size.l;
    if (il_ar >= lut->ar_size.l)
        il_ar = lut->ar_size.l - 1;

    /* Loop through the samples in this line */
    for (is = 0; is < nsamp; is++) {
        loc.s = is;
        cld_col = is / cld_diags->cellwidth;
        is_ar = is / lut->ar_region_size.s;
        if (is_ar >= lut->ar_size.s)
            is_ar = lut->ar_size.s - 1;

        is_fill = false;
        if (thermal_band) {
            if (b6_line[is] == lut->b6_in_fill) {
                is_fill = true;
                ddv_line[is] = 0x08;
            }
        }

        if ((qa_line[is] & 0x01) == 0x01) {
            is_fill = true;
            ddv_line[is] = 0x08;
        }

        if (! is_fill) {
            ddv_line[is] &= 0x44; /* reset all bits except cloud shadow and
                                     adjacent cloud */ 

            if (((qa_line[is] & 0x08) == 0x08) ||
                ((lut->meta.inst == INST_TM) &&
                (line_in[2][is] >= thresh_TM))) {
                if (thermal_band) {
                    t6 = b6_line[is] * lut->b6_scale_factor + lut->b6_add_offset;

                    /* Interpolate the cloud diagnostics for current pixel */
                    interpol_clddiags_1pixel (cld_diags, il, is, tmpflt_arr);
                    temp_b6_clear = tmpflt_arr[0];
                    atemp_ancillary = tmpflt_arr[1];
                    if (temp_b6_clear < 0.) {
                        temp_thshld1 = atemp_ancillary - 20.;
                        temp_thshld2 = atemp_ancillary - 20.;
                    }
                    else {
                        if (cld_diags->std_t6_clear[cld_row][cld_col] > 0.) {
                            temp_thshld1 = temp_b6_clear -
                               (cld_diags->std_t6_clear[cld_row][cld_col] + 4.);
                            temp_thshld2 = temp_b6_clear -
                               cld_diags->std_t6_clear[cld_row][cld_col];
                        }
                        else {
                            temp_thshld1 = temp_b6_clear - 4.;
                            temp_thshld2 = temp_b6_clear - 2.;
                        }
                    }

                    if ((((qa_line[is] & 0x20) == 0x20) ||
                         ((lut->meta.inst == INST_TM) &&
                          (line_in[4][is] >= thresh_TM))) && (t6 < temp_thshld1)) {
                        /* saturated band 5 and t6 < threshold => cloudy */
                        ddv_line[is] &= 0xbf; /* reset shadow bit */
                        ddv_line[is] &= 0xfb; /* reset adjacent cloud bit */
                        ddv_line[is] |= 0x20; /* set cloud bit */
                    }
                    else if ((line_in[4][is] < band5_thresh) &&
                             (t6 < temp_snow_thshld)) { /* snow */
                        ddv_line[is] |= 0x80;
                    }
                    else { /* assume cloudy */
                        ddv_line[is] &= 0xbf; /* reset shadow bit */
                        ddv_line[is] &= 0xfb; /* reset adjacent cloud bit */
                        ddv_line[is] |= 0x20; /* set cloud bit */
                    }
                }
            }
            else {
                /* Interpolate the atmospheric conditions for current pixel */
                SrInterpAtmCoef (lut, &loc, &atmos_coef, &interpol_atmos_coef);

                /* Compute the reflectance for each band using the interpolated
                   atmospheric coefficients */
                rho1 = line_in[0][is] * lut->scale_factor + lut->add_offset;
                rho1 = (rho1/interpol_atmos_coef.tgOG[0][0] -
                    interpol_atmos_coef.rho_ra[0][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[0][0] *
                    interpol_atmos_coef.td_ra[0][0] *
                    interpol_atmos_coef.tu_ra[0][0]);
                rho1 /= tmpflt;
                rho1 /= (1. + interpol_atmos_coef.S_ra[0][0] * rho1);

                rho2 = line_in[1][is] * lut->scale_factor + lut->add_offset;
                rho2 = (rho2/interpol_atmos_coef.tgOG[1][0] -
                    interpol_atmos_coef.rho_ra[1][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[1][0] *
                    interpol_atmos_coef.td_ra[1][0] *
                    interpol_atmos_coef.tu_ra[1][0]);
                rho2 /= tmpflt;
                rho2 /= (1. + interpol_atmos_coef.S_ra[1][0] * rho2);

                rho3 = line_in[2][is] * lut->scale_factor + lut->add_offset;
                rho3 = (rho3 / interpol_atmos_coef.tgOG[2][0] -
                    interpol_atmos_coef.rho_ra[2][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[2][0] *
                    interpol_atmos_coef.td_ra[2][0] *
                    interpol_atmos_coef.tu_ra[2][0]);
                rho3 /= tmpflt;
                rho3 /= (1. + interpol_atmos_coef.S_ra[2][0] * rho3);

                rho4 = line_in[3][is] * lut->scale_factor + lut->add_offset;
                rho4 = (rho4 / interpol_atmos_coef.tgOG[3][0] -
                    interpol_atmos_coef.rho_ra[3][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[3][0] *
                    interpol_atmos_coef.td_ra[3][0] *
                    interpol_atmos_coef.tu_ra[3][0]);
                rho4 /= tmpflt;
                rho4 /= (1. + interpol_atmos_coef.S_ra[3][0] * rho4);

                rho5 = line_in[4][is] * lut->scale_factor + lut->add_offset;
                rho5 = (rho5 / interpol_atmos_coef.tgOG[4][0] -
                    interpol_atmos_coef.rho_ra[4][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[4][0] *
                    interpol_atmos_coef.td_ra[4][0] *
                    interpol_atmos_coef.tu_ra[4][0]);
                rho5 /= tmpflt;
                rho5 /= (1. + interpol_atmos_coef.S_ra[4][0] * rho5);

                rho7 = line_in[5][is] * lut->scale_factor + lut->add_offset;
                rho7 = (rho7 / interpol_atmos_coef.tgOG[5][0] -
                    interpol_atmos_coef.rho_ra[5][0]);
                tmpflt = (interpol_atmos_coef.tgH2O[5][0] *
                    interpol_atmos_coef.td_ra[5][0] *
                    interpol_atmos_coef.tu_ra[5][0]);
                rho7 /= tmpflt;
                rho7 /= (1. + interpol_atmos_coef.S_ra[5][0] * rho7);

                /* Get the temperature */
                if (thermal_band)
                    t6 = b6_line[is] * lut->b6_scale_factor + lut->b6_add_offset;

                /* Interpolate the cloud diagnostics for the current pixel */
                interpol_clddiags_1pixel (cld_diags, il, is, tmpflt_arr);
                temp_b6_clear = tmpflt_arr[0];
                atemp_ancillary = tmpflt_arr[1];

                if (temp_b6_clear < 0.) {
                    temp_thshld1 = atemp_ancillary - 20.;
                    temp_thshld2 = atemp_ancillary - 20.;
                }
                else {
                    if (cld_diags->std_t6_clear[cld_row][cld_col] > 0.) {
                        temp_thshld1 = temp_b6_clear -
                            (cld_diags->std_t6_clear[cld_row][cld_col] + 4.);
                        temp_thshld2 = temp_b6_clear -
                            cld_diags->std_t6_clear[cld_row][cld_col];
                    }
                    else {
                        temp_thshld1 = temp_b6_clear - 4.;
                        temp_thshld2 = temp_b6_clear - 2.;
                    }
                }

                if (thermal_band) {
                    /* Compute cloud coefficients */
                    vra = rho1 - rho3 * 0.5;

                    C1 = (int)(vra > VRA_THRESHOLD);
                    C2 = (t6 < temp_thshld1);
   
                    tmpflt = rho4 / rho3;

                    C4 = (rho7 > 0.03);
                    C5 = (t6 < temp_thshld2) && C1;
                }

                /**
                Water test :
                ndvi < 0 => water
                ((0<ndvi<0.1) or (b4<5%)) and b5 < 0.01 => turbid water
                **/
                if ((rho4 + rho3) != 0)
                    ndvi = (rho4 - rho3) / (rho4 + rho3);
                else
                    ndvi = 0.01;
                water = (ndvi < 0) || ((((ndvi > 0) && (ndvi < 0.1)) ||
                    (rho4 < 0.05)) && (rho5 < 0.02));

                if (thermal_band) {
                    if (!water) { /* if not water */
                        ddv_line[is] |= 0x10;
                        if ((C2 || C5) && C4) { /* cloudy */
                            ddv_line[is] &= 0xbf; /* reset shadow bit */
                            ddv_line[is] &= 0xfb; /* reset adjacent cloud bit */
                            ddv_line[is] |= 0x20; /* set cloud bit */
                        }
                        else { /* clear */
                            ddv_line[is] &= 0xdf;
                            ndsi = (rho2 - rho5) / (rho2 + rho5);
                            if ((ndsi > 0.3) && (t6 < temp_snow_thshld) &&
                                (rho4 > 0.2))
                                ddv_line[is] |= 0x80;
                        }
                    }
                    else 
                        ddv_line[is] &= 0xef; 
                }
                else { /* no thermal band - cannot run cloud mask */
                    ddv_line[is] &= 0xdf; /* assume clear */
                    if (!water) { /* if not water */
                        ddv_line[is] |= 0x10;
                    }
                    else {
                        ddv_line[is] &= 0xef; 
                    }
                }
            }  /* end else saturated band 3 */
        }  /* if ! is_fill */ 
    }  /* end for is */

    free_mem_atmos_coeff(&interpol_atmos_coef);
    return true;
}


bool dilate_cloud_mask
(
    Lut_t *lut,
    int nsamp,
    char ***cloud_buf,
    int dilate_dist
)
{
    int il,is,il_adj,is_adj,buf_ind;
    int k;

    for (il = 0; il < lut->ar_region_size.l; il++) {
        for (is = 0; is < nsamp; is++) {
            if (cloud_buf[1][il][is] & 0x20) { /* if cloudy dilate */
                for (k = (il-dilate_dist); k < (il+dilate_dist); k++) {
                    il_adj = k;
                    buf_ind = 1;
                    if (k < 0) {
                        buf_ind--;
                        il_adj += lut->ar_region_size.l; 
                    }
                    if (k >= lut->ar_region_size.l) {
                        buf_ind++;
                        il_adj -= lut->ar_region_size.l; 
                    }
                    if ((il_adj>=0)&&(il_adj<lut->ar_region_size.l)) {
                        for (is_adj = (is-dilate_dist);
                             is_adj < (is+dilate_dist);
                             is_adj++) {
                            if ((is_adj>=0)&&(is_adj<nsamp)) {
                                /* if not cloudy */
                                if (!(cloud_buf[buf_ind][il_adj][is_adj] &
                                    0x20)) {
                                    /* reset adjacent cloud bit */
                                    cloud_buf[buf_ind][il_adj][is_adj] &= 0xfb;
                                    /* reset shadow bit */
                                    cloud_buf[buf_ind][il_adj][is_adj] &= 0xbf;
                                    /* set adjacent cloud bit */
                                    cloud_buf[buf_ind][il_adj][is_adj] |= 0x04;
                                }
                            }
                        }  /* for is_adj */
                    }  /* if il_adj */
                }  /* for k */
            } /* if cloudy */
        }  /* for is */
    }  /* for il */
    return true;
}


void cast_cloud_shadow
(
    Lut_t *lut,
    int nsamp,
    int il_start,
    uint16_t ***line_in,
    uint16_t **b6_line,
    cld_diags_t *cld_diags,
    char ***cloud_buf,
    Ar_gridcell_t *ar_gridcell,
    float pixel_size,
    float adjust_north
)
{
    int il,is,il_ar,is_ar,shd_buf_ind;
    float t6,temp_b6_clear,atemp_ancillary,tmpflt_arr[10];
    float conv_factor,cld_height,ts,tv,fs,fv,dx,dy;
    int shd_x,shd_y;

/***
    Cloud Shadow
***/
    il_ar = il_start / lut->ar_region_size.l;
    if (il_ar >= lut->ar_size.l)
        il_ar = lut->ar_size.l - 1;
    for (il = 0; il <lut->ar_region_size.l; il++) {
        for (is = 0; is < nsamp; is++) {
            is_ar = is / lut->ar_region_size.s;
            if (is_ar >= lut->ar_size.s)
                is_ar = lut->ar_size.s - 1;

            /* Get the thermal info */
            t6 = b6_line[il][is] * lut->b6_scale_factor + lut->b6_add_offset;

            /* Interpolate the cloud diagnostics for this pixel */
            interpol_clddiags_1pixel (cld_diags, il+il_start, is, tmpflt_arr);
            temp_b6_clear = tmpflt_arr[0];
            atemp_ancillary = tmpflt_arr[2];

            if (cloud_buf[1][il][is] & 0x20) { /* if cloudy cast shadow */
                conv_factor = 6.;
                while (conv_factor <= 6.) {
                    /* Determine the cloud height */
                    if (temp_b6_clear > 0)
                        cld_height = (temp_b6_clear - t6) / conv_factor;
                    else
                        cld_height = (atemp_ancillary - t6) / conv_factor;

                    /* If the cloud height is greater than 0, then determine
                       the shadow */
                    if (cld_height > 0.) {
                        ts = ar_gridcell->sun_zen[il_ar*lut->ar_size.s+is_ar]
                            / DEG;
                        fs = (ar_gridcell->rel_az[il_ar*lut->ar_size.s+is_ar]
                            - adjust_north) / DEG;
                        tv = ar_gridcell->view_zen[il_ar*lut->ar_size.s+is_ar]
                            / DEG;
                        fv = 0.;

/*                        dy = sin(2.*M_PI - fv) * tan(tv) * cld_height;
                        dx = cos(2.*M_PI - fv) * tan(tv) * cld_height;
                        dy=0;
                        dx=0;

                        shd_x = is - dx * 1000. / pixel_size;
                        shd_y = il + dy * 1000. / pixel_size;
*/
                        dy = cos(fs) * tan(ts) * cld_height;
                        dx = sin(fs) * tan(ts) * cld_height;
                        shd_x = is - dx * 1000. / pixel_size;
                        shd_y = il + dy * 1000. / pixel_size;

                        if ((shd_x >= 0) && (shd_x < nsamp)) {
                            shd_buf_ind = 1;
                            if (shd_y < 0) {
                                shd_buf_ind--;
                                shd_y += lut->ar_region_size.l;
                            }
                            if (shd_y >= lut->ar_region_size.l) {
                                shd_buf_ind++;
                                shd_y -= lut->ar_region_size.l;
                            }
                            /* Mask as cloud shadow */
                            if (shd_y >= 0 && shd_y < lut->ar_region_size.l) {
                                /* if not cloud, adjacent cloud or cloud
                                   shadow */
                                if (!((cloud_buf[shd_buf_ind][shd_y][shd_x] &
                                       0x20) ||
                                      (cloud_buf[shd_buf_ind][shd_y][shd_x] &
                                       0x04) ||
                                      (cloud_buf[shd_buf_ind][shd_y][shd_x] &
                                       0x40)))
                                    /* set cloud shadow bit */
                                   cloud_buf[shd_buf_ind][shd_y][shd_x] |= 0x40;
                            }
                        }
                    } /* if cld_height > 0 */

                    conv_factor += 1.;
                } /* while conv_fact <= 6. */
            }
        }
    }

    return;
}

bool dilate_shadow_mask
(
    Lut_t *lut,          /* I: lookup table */
    int nsamp,           /* I: number of samples in the current line */
    char ***cloud_buf,   /* I: I/O: cloud buffer */
    int dilate_dist      /* I: size of dilation window */
)
{
    int il,is,il_adj,is_adj,buf_ind;
    int k;
    char *fill_mask = NULL;

    if ((fill_mask = calloc(lut->ar_region_size.l*nsamp, sizeof(char))) == NULL)
        return false;

    for (il = 0; il < lut->ar_region_size.l; il++) {
        for (is = 0; is < nsamp; is++) {
            if ((cloud_buf[0][il][is] & 0x40) && (!fill_mask[il*nsamp+is])) {
                /* if cloud shadow dilate */
                for (k = (il-dilate_dist); k <= (il+dilate_dist); k++) {
                    il_adj = k;
                    buf_ind = 0;
                    if (k >= lut->ar_region_size.l) {
                        buf_ind++;
                        il_adj -= lut->ar_region_size.l; 
                    }

                    if (k >= 0) {
                        if ((il_adj >= 0) && (il_adj < lut->ar_region_size.l)) {
                            for (is_adj = (is-dilate_dist);
                                 is_adj <= (is+dilate_dist); is_adj++) {
                                if ((is_adj >= 0) && (is_adj < nsamp)) {
                                    /* if not cloud, adjacent cloud or cloud
                                       shadow */
                                    if (!((cloud_buf[buf_ind][il_adj][is_adj] &
                                           0x20) ||
                                          (cloud_buf[buf_ind][il_adj][is_adj] &
                                           0x04) ||
                                          (cloud_buf[buf_ind][il_adj][is_adj] &
                                           0x40))) {
                                        /* set adjacent cloud shadow bit */
                                        cloud_buf[buf_ind][il_adj][is_adj] |=
                                            0x40;
                                        fill_mask[il_adj*nsamp+is_adj]=1;
                                    }
                                }
                            }
                        }
                    }
                }
            } /* if cloud shadow */
        }  /* for is */
    }  /* for il */

    free (fill_mask);
    return true;
}

int allocate_cld_diags
(
    cld_diags_t *cld_diags,   /* O: cloud diagnostics */
    int cell_height,          /* I: cell height of cloud diagnostics */
    int cell_width,           /* I: cell width of cloud diagnostics */
    int scene_height,         /* I: number of lines in the scene */
    int scene_width           /* I: number of samples in the scene */
)
{
    int i;
    
    /* Set the size of the diagnostic cells and the number of rows/cols in
       the diagnostics */
    cld_diags->cellheight = cell_height;
    cld_diags->cellwidth = cell_width;
    cld_diags->nbrows = (scene_height - 1) / cell_height + 1;
    cld_diags->nbcols = (scene_width - 1) / cell_width + 1;

    if ((cld_diags->avg_t6_clear = malloc (cld_diags->nbrows*sizeof(float *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++)
        if ((cld_diags->avg_t6_clear[i] = calloc(cld_diags->nbcols,
            sizeof(float)))==NULL)
            return -1;

    if ((cld_diags->std_t6_clear = malloc (cld_diags->nbrows*sizeof(float *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++)
        if ((cld_diags->std_t6_clear[i] = calloc(cld_diags->nbcols,
            sizeof(float)))==NULL)
            return -1;

    if ((cld_diags->avg_b7_clear = malloc (cld_diags->nbrows*sizeof(float *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->avg_b7_clear[i] = calloc(cld_diags->nbcols,
            sizeof(float)))==NULL)
            return -1;

    if ((cld_diags->std_b7_clear = malloc (cld_diags->nbrows*sizeof(float *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->std_b7_clear[i] = calloc(cld_diags->nbcols,
            sizeof(float)))==NULL)
            return -1;

    if ((cld_diags->airtemp_2m = malloc (cld_diags->nbrows*sizeof(float *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->airtemp_2m[i] = calloc(cld_diags->nbcols,
            sizeof(float)))==NULL)
            return -1;

    if ((cld_diags->nb_t6_clear = malloc (cld_diags->nbrows*sizeof(int *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->nb_t6_clear[i] = calloc(cld_diags->nbcols,
            sizeof(int)))==NULL)
            return -1;

    return 0;
}

void free_cld_diags(cld_diags_t *cld_diags) {
    
    int i;
    
    for (i = 0; i < cld_diags->nbrows; i++) {
        free(cld_diags->avg_t6_clear[i]);
        free(cld_diags->std_t6_clear[i]);
        free(cld_diags->avg_b7_clear[i]);
        free(cld_diags->std_b7_clear[i]);
        free(cld_diags->airtemp_2m[i]);
        free(cld_diags->nb_t6_clear[i]);
    }

    free(cld_diags->avg_t6_clear);
    free(cld_diags->std_t6_clear);
    free(cld_diags->avg_b7_clear);
    free(cld_diags->std_b7_clear);
    free(cld_diags->airtemp_2m);
    free(cld_diags->nb_t6_clear);
}

void fill_cld_diags(cld_diags_t *cld_diags) {
/*
!Description: fill in missing values in the T6Clear grid based
on existing values (spatial interpolation). Missing values have been previously
set to -9999. A filling can be distincted from a original value by looking at
the standard deviation of the optical depth which is set to -9999 for a filling.
!END****************************************************************************
*/
    int i,j,k,l,pass,count;
    float lastt6=0.0, dist, sumt6=0.0, sumdist=0.0, lastb7=0.0, sumb7=0.0;
    char missing_flag[300][300];
    int min_nb_values,n,max_distance;
   
    count=0;
    for (i=0;i<cld_diags->nbrows;i++) {
        for (j=0;j<cld_diags->nbcols;j++) {
            missing_flag[i][j]=1;
            if (cld_diags->avg_t6_clear[i][j]!=-9999.)  {
                count++;
                lastt6=cld_diags->avg_t6_clear[i][j];
                lastb7=cld_diags->avg_b7_clear[i][j];
                missing_flag[i][j]=0;
            }
        }
    }

    if (count==0)
        return;

    else if (count==1) {
        for (i=0;i<cld_diags->nbrows;i++)
            for (j=0;j<cld_diags->nbcols;j++) {
                cld_diags->avg_t6_clear[i][j]=lastt6;
                cld_diags->avg_b7_clear[i][j]=lastb7;
            }

        return;
    }
    
    /* Loop through the lines and samples in the cloud diagnostics */
    for (i=0;i<cld_diags->nbrows;i++) {
        for (j=0;j<cld_diags->nbcols;j++) {
            /**
            Look for at least 3 neighboring valid values within 4 GPs
            **/
            min_nb_values=3;
            max_distance=4;    
            pass=0;
            while ((cld_diags->avg_t6_clear[i][j] == -9999.) &&
                   (pass<max_distance)) {
                pass++;
                sumdist=0.;
                sumt6=0.;
                sumb7=0.;
                n=0;
                for (k=i-pass;k<=(i+pass);k++) {
                    if ((k>=0)&&(k<cld_diags->nbrows)) {
                        for (l=j-pass;l<=(j+pass);l++) {
                            if ((l>=0)&&(l<cld_diags->nbcols)) {
                                if (!missing_flag[k][l]) {
                                    dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                    sumdist += dist;
                                    sumt6 += (dist *
                                        cld_diags->avg_t6_clear[k][l]);
                                    sumb7 += (dist *
                                        cld_diags->avg_b7_clear[k][l]);
                                    n++;
                                }
                            }
                        }
                    }
                }

                if ((n>=min_nb_values)&&(sumdist!=0.)) {
                    cld_diags->avg_t6_clear[i][j] = sumt6 / sumdist;
                    cld_diags->avg_b7_clear[i][j] = sumb7 / sumdist;
                }
            }  /* end while */

            /**
            Look for at least 2 neighboring valid values within 6 GPs
            **/
            min_nb_values=2;
            max_distance=6;    
            pass=0;
            while ((cld_diags->avg_t6_clear[i][j] == -9999.) &&
                   (pass<max_distance)) {
                pass++;
                sumdist=0.;
                sumt6=0.;
                sumb7=0.;
                n=0;
                for (k=i-pass;k<=(i+pass);k++) {
                    if ((k>=0)&&(k<cld_diags->nbrows)) {
                        for (l=j-pass;l<=(j+pass);l++) {
                            if ((l>=0)&&(l<cld_diags->nbcols)) {
                                if (!missing_flag[k][l]) {
                                    dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                    sumdist += dist;
                                    sumt6 += (dist *
                                        cld_diags->avg_t6_clear[k][l]);
                                    sumb7 += (dist *
                                        cld_diags->avg_b7_clear[k][l]);
                                    n++;
                                }
                            }
                        }
                    }
                }

                if ((n>=min_nb_values)&&(sumdist!=0.)) {
                    cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
                    cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
                }
            }  /* end while */

            /**
            Look for at least 1 neighboring valid values within 10 GPs
            **/
            min_nb_values=1;
            max_distance=10;    
            pass=0;
            while ((cld_diags->avg_t6_clear[i][j] == -9999.) &&
                   (pass<max_distance)) {
                pass++;
                sumdist=0.;
                sumt6=0.;
                sumb7=0.;
                n=0;
                for (k=i-pass;k<=(i+pass);k++) {
                    if ((k>=0)&&(k<cld_diags->nbrows)) {
                        for (l=j-pass;l<=(j+pass);l++) {
                            if ((l>=0)&&(l<cld_diags->nbcols)) {
                                if (!missing_flag[k][l]) {
                                    dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                    sumdist += dist;
                                    sumt6 += (dist *
                                        cld_diags->avg_t6_clear[k][l]);
                                    sumb7 += (dist *
                                        cld_diags->avg_b7_clear[k][l]);
                                    n++;
                                }
                            }
                        }
                    }
                }

                if ((n>=min_nb_values)&&(sumdist!=0.)) {
                    cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
                    cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
                }
            }  /* end while */
        }  /* for j */
    }  /* for i */
}

void interpol_clddiags_1pixel
(
    cld_diags_t *cld_diags,  /* I: cloud diagnostics */
    int img_line,            /* I: current line in image */
    int img_sample,          /* I: current sample in image */
    float *inter_value       /* O: interpolated cloud diagnostic value */
)
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

    inter_value[0] => t6_clear
    inter_value[1] => airtemp_2m

    Updated by Gail Schmidt, USGS EROS, on 10/20/2014
    We want the airtemp_2m to be calculated regardless of whether the band6
        temp is available (i.e. there were clear pixels to compute the average
        thermal temp).  Many of the users of these interpolated values use the
        airtemp_2m as the default if the band6 clear temps are not valid.
 */

{
    typedef struct {
      int l;                /* line number */
      int s;                /* sample number */
    } Img_coord_int_t;
    
    Img_coord_int_t p[4];
    int i, n, n_anc;
    float dl, ds, w;
    float sum[10], sum_w, sum_anc_w;

    int cell_half_height, cell_half_width;

    for (i = 0; i < 3; i++) 
        inter_value[i] = -9999.;

    cell_half_height = (cld_diags->cellheight + 1) >> 1;  /* divide by 2 */
    cell_half_width = (cld_diags->cellwidth + 1) >> 1;  /* divide by 2 */

    p[0].l = (img_line - cell_half_height) / cld_diags->cellheight;
    if (p[0].l < 0)
        p[0].l = 0;
    p[2].l = p[0].l + 1;
    if (p[2].l >= cld_diags->nbrows) {
        p[2].l = cld_diags->nbrows - 1;
            if (p[0].l > 0)
                p[0].l--;
    }
        
    p[1].l = p[0].l;
    p[3].l = p[2].l;

    p[0].s = (img_sample - cell_half_width) / cld_diags->cellwidth;
    if (p[0].s < 0)
        p[0].s = 0;
    p[1].s = p[0].s + 1;

    if (p[1].s >= cld_diags->nbcols) {
        p[1].s = cld_diags->nbcols - 1;
        if (p[0].s > 0)
            p[0].s--;
    }

    p[2].s = p[0].s;
    p[3].s = p[1].s;

    /* Initialize the counting and sum variables */
    n = 0;
    n_anc = 0;
    sum_w = 0.0;
    sum_anc_w = 0.0;
    for (i = 0; i < 2; i++)
        sum[i]=0.;

    /* Loop through the four points to be used in the interpolation */
    for (i = 0; i < 4; i++) {
        /* If the points are valid */
        if (p[i].l != -1 && p[i].s != -1) {
            dl = fabs(img_line - cell_half_height) -
                (p[i].l * cld_diags->cellheight);
            dl = fabs(dl) / cld_diags->cellheight; 
            ds = fabs(img_sample - cell_half_width) -
                (p[i].s * cld_diags->cellwidth);
            ds = fabs(ds) / cld_diags->cellwidth; 
            w = (1.0 - dl) * (1.0 - ds);

            if (cld_diags->avg_t6_clear[p[i].l][p[i].s] != -9999.) {
                n++;
                sum_w += w;
                sum[0] += (cld_diags->avg_t6_clear[p[i].l][p[i].s] * w);
            }

            if (cld_diags->airtemp_2m[p[i].l][p[i].s] != -9999) {
                n_anc++;
                sum_anc_w += w;
                sum[1] += (cld_diags->airtemp_2m[p[i].l][p[i].s] * w);
            }
        }  /* end if points are valid */
    }  /* end for i */

    if ((n > 0) && (sum_w > 0)) {
        inter_value[0] = sum[0] / sum_w;
    }

    if ((n_anc > 0) && (sum_anc_w > 0)) {
        inter_value[1] = sum[1] / sum_anc_w;
    }

    return;
}
