#ifndef LNDSR_H
#define LNDSR_H

#define bounded(A,B,C) (A>B?(A<C?A:C):B)
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))

#include <stdio.h>
#include "hdf.h"
#include "mfhdf.h"
#include "mystring.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "espa_geoloc.h"

/* Extra bands - atmos_opacity, cloud_QA */
#define NBAND_SR_EXTRA (2)
#define NBAND_REFL_MAX (6)
#define NBAND_PRWV_MAX (3)
#define NBAND_SR_MAX (NBAND_REFL_MAX + NBAND_SR_EXTRA)

typedef enum {
  ATMOS_OPACITY = 0,
  CLOUD,
  FILL,
  DDV,
  CLOUD_SHADOW,
  SNOW,
  LAND_WATER,
  ADJ_CLOUD
} Extra_Band_t;

typedef enum {
  QA_OFF = 0,
  QA_ON = 255
} QA_t;

/* Bit representations for the bit-packed QA band */
typedef enum {
  DDV_BIT = 0,
  CLOUD_BIT,
  CLOUD_SHADOW_BIT,
  ADJ_CLOUD_BIT,
  SNOW_BIT,
  LAND_WATER_BIT,
} QA_Band_Bit_t;

/* Satellite type definition */

typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_1 = 0, 
  SAT_LANDSAT_2, 
  SAT_LANDSAT_3, 
  SAT_LANDSAT_4, 
  SAT_LANDSAT_5, 
  SAT_LANDSAT_7, 
  SAT_MAX
} Sat_t;

extern const Key_string_t Sat_string[SAT_MAX];

/* Instrument type definition */

typedef enum {
  INST_NULL = -1,
  INST_MSS = 0, 
  INST_TM,
  INST_ETM, 
  INST_MAX
} Inst_t;

extern const Key_string_t Inst_string[INST_MAX];

/* World Reference System (WRS) type definition */

typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
  WRS_MAX
} Wrs_t;

extern const Key_string_t Wrs_string[WRS_MAX];

/* Ozone source */

typedef enum {
  OZSRC_NULL = -1,
  OZSRC_NIMBUS7 = 0, 
  OZSRC_METEOR3, 
  OZSRC_EARTHPROBE, 
  OZSRC_OMI, 
  OZSRC_NIMBUS7_FILL, 
  OZSRC_METEOR3_FILL, 
  OZSRC_EARTHPROBE_FILL, 
  OZSRC_FILL, 
  OZSRC_MAX
} Ozsrc_t;

extern const Key_string_t Ozsrc_string[OZSRC_MAX];

typedef struct {
  int nbrows,nbcols;
  float *lat,*lon;
  float *sun_zen,*view_zen,*rel_az;
  float *wv,*spres,*ozone,*spres_dem;
  float *line_lat,*line_lon,*line_sun_zen,*line_view_zen,*line_rel_az;
  float *line_wv,*line_spres,*line_ozone,*line_spres_dem;
} Ar_gridcell_t;

typedef struct {
int *computed;
float *tgOG[7],*tgH2O[7],*td_ra[7],*tu_ra[7],*rho_mol[7],*rho_ra[7],*td_da[7],*tu_da[7],*S_ra[7];
float *td_r[7],*tu_r[7],*S_r[7],*rho_r[7];
} atmos_t;

int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef);
int free_mem_atmos_coeff(atmos_t *atmos_coef);

#endif
