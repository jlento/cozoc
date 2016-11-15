#ifndef CONTEXT_H
#define CONTEXT_H

#include "petscksp.h"

enum {
        OMEGA_QUASI_GEOSTROPHIC = (1 << 0),
        OMEGA_GENERALIZED       = (1 << 1),
        HEIGHT_TENDENCY         = (1 << 2)
};

typedef struct {
        DM           da,da2;
        DM           daxy;
        KSP          ksp;
        int          ncid;      /* NetCDF file id */
        int          eqns;      /* The equations to be calculated */
        size_t       mx,my,mz,mt; /* Global grid sizes */
        PetscScalar  hx,hy,hz;  /* Grid spacings */
        //        Vec          p;         /* Pressure, p[mz]*/
        Vec          f;         /* Coriolis parameter, f[my]) */
        PetscInt     time;      /* Current timestep */
        Vec          psfc;      /* Surface pressure */
        Vec          T;         /* Temperature */
        Vec          Z;         /* Geopotential height */
        PetscScalar *Pressure;
        PetscScalar *Coriolis_parameter;
        Vec          Sigma_parameter;
        Vec          Vorticity;
        Vec          Horizontal_wind;

} tContext,*Context;

extern PetscErrorCode context_create(const char *fname,int *eqns,Context *ctx);

extern PetscErrorCode context_destroy(Context *ctx);

extern PetscErrorCode context_update(const int time,Context ctx);

#endif /* CONTEXT_H */
