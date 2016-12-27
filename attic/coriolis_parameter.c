#include "context.h"
#include "coriolis_parameter.h"
#include "petscdmda.h"
#include "utils.h"


int coriolis_parameter_add (Vec x, Context ctx) {

    PetscScalar*   f = ctx->Coriolis_parameter;
    DM             da;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar*** xa;

    VecGetDM (x, &da);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);
    DMDAVecGetArray (da, x, &xa);

    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                xa[k][j][i] += f[j]; } } }

    DMDAVecRestoreArray (da, x, &xa);
    return (0); }
