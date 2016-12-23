#include "context.h"
#include "coriolis_parameter.h"
#include "petscdmda.h"
#include "utils.h"


extern PetscErrorCode coriolis_parameter_add (Vec x, Context ctx) {
    PetscScalar*   f = ctx->Coriolis_parameter;
    DM             da;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar*** xa;
    PetscErrorCode ierr;

    ierr = VecGetDM (x, &da);
    ierr = DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);
    ierr = DMDAVecGetArray (da, x, &xa);

    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                xa[k][j][i] += f[j]; } } }

    ierr = DMDAVecRestoreArray (da, x, &xa);
    return (0); }
