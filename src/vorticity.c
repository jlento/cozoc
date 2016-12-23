#include "constants.h"
#include "derivatives.h"
#include "loops.h"
#include "petscdmda.h"
#include "vorticity.h"

extern PetscErrorCode vorticity (Context ctx) {

    DM              da       = ctx->da;
    DM              da2      = ctx->da2;
    Vec             Vvec     = ctx->Horizontal_wind;
    Vec             sigmavec = ctx->Sigma_parameter;
    Vec             zetavec  = ctx->Vorticity;
    PetscScalar     hz       = ctx->hz;
    size_t          mz       = ctx->mz;
    PetscScalar*    f        = ctx->Coriolis_parameter;
    PetscScalar***  sigma;
    PetscScalar**** V;
    PetscScalar***  zeta;
    PetscScalar     zetap, dudp, dvdp, dVdp2;

    horizontal_rotor (Vvec, zetavec, ctx);

    DMDAVecGetArray (da, zetavec, &zeta);
    DMDAVecGetArrayRead (da, sigmavec, &sigma);
    DMDAVecGetArrayDOFRead (da2, Vvec, &V);

    LOOP_KJI (da, if (f[j] > 1.0e-7) {
    int         k0;
    int         k1;
    PetscScalar w;

    if (k == 0) {
            k1 = 1;
            k0 = 0;
            w  = 1.0 / (hz * hz); }
        else if (k == mz - 1) {
            k1 = mz - 1;
            k0 = mz - 2;
            w  = 1.0 / (hz * hz); }
        else {
            k1 = k + 1;
            k0 = k - 1;
            w  = 0.25 / (hz * hz); }

        dudp  = w * (V[k1][j][i][0] - V[k0][j][i][0]);
        dvdp  = w * (V[k1][j][i][1] - V[k0][j][i][1]);
        dVdp2 = dudp * dudp + dvdp * dvdp;
        zetap = etamin + f[j] / (4.0 * sigma[k][j][i]) * dVdp2 - f[j];

        if (zetap > zeta[k][j][i])
            zeta[k][j][i] = zetap; });
    DMDAVecRestoreArray (da, zetavec, &zeta);
    DMDAVecRestoreArrayRead (da, sigmavec, &sigma);
    DMDAVecRestoreArrayDOFRead (da2, Vvec, &V);

    return (0); }
