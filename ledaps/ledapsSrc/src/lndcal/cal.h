#ifndef CAL_H
#define CAL_H

#include <stdbool.h>
#include "lndcal.h"
#include "lut.h"
#include "input.h"
static const int SATU_VAL[7]={255,255,255,255,255,255,255};
static const int SATU_VAL6= 254;

typedef struct {
  bool first[NBAND_REFL_MAX];
  unsigned char idn_min[NBAND_REFL_MAX];
  unsigned char idn_max[NBAND_REFL_MAX];
  float rad_min[NBAND_REFL_MAX];
  float rad_max[NBAND_REFL_MAX];
  float ref_min[NBAND_REFL_MAX];
  float ref_max[NBAND_REFL_MAX];
  int iref_min[NBAND_REFL_MAX];
  int iref_max[NBAND_REFL_MAX];
} Cal_stats_t;

typedef struct {
  bool first;
  unsigned char idn_min;
  unsigned char idn_max;
  float rad_min;
  float rad_max;
  float temp_min;
  float temp_max;
  int itemp_min;
  int itemp_max;
} Cal_stats6_t;

bool Cal(Param_t *param, Lut_t *lut, int iband, Input_t *input,
         unsigned char *line_in, int16 *line_in_sun_zen, uint16_t *line_out,
         unsigned char *line_out_qa, Cal_stats_t *cal_stats, int iy);

bool Cal6(Lut_t *lut, Input_t *input, unsigned char *line_in,
  uint16_t *line_out, unsigned char *line_out_qa, Cal_stats6_t *cal_stats, int iy);

#endif
