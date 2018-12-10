#ifndef AR_H
#define AR_H

#include <stdbool.h>
#include "lndsr.h"
#include "lut.h"
#include "sixs_runs.h"

typedef struct {
  bool first;
  int ar_min;
  int ar_max;
  long nfill;
} Ar_stats_t;

typedef struct
{
    short b[3];   /* collect_band */
    short b7;     /* collect_band7 */
} collect_bands_t;

bool Ar(int il_ar, Lut_t *lut, Img_coord_int_t *size_in, int16 ***line_in,
        char **ddv_line, atmos_t *atmos_coef_ar, collect_bands_t *cbands,
        int **line_ar, Ar_stats_t *ar_stats, Ar_gridcell_t *ar_gridcell,
        sixs_tables_t *sixs_tables);

int ArInterp(Lut_t *lut, Img_coord_int_t *loc, int ***line_ar, int *inter_aot);
int Fill_Ar_Gaps(Lut_t *lut, int ***line_ar, int ib);

#endif
