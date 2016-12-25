#include "constants.h"
#include "derivatives.h"
#include "loops.h"
#include "petscdmda.h"
#include "sigma_parameter.h"

int sigma_parameter (Context ctx) {
    DM             da       = ctx->da;
    PetscScalar    hz       = ctx->hz;
    PetscScalar*   p        = ctx->Pressure;
    Vec            Tvec     = ctx->Temperature;
    Vec            sigmavec = ctx->Sigma_parameter;
    const double   R        = Specific_gas_constant_of_dry_air;
    const double   c_p      = Specific_heat_of_dry_air;
    PetscScalar*** T;
    PetscScalar*** sigma;

    /* first sigma holds auxiliary variable dT/dp, temporarily */
    pder (hz, Tvec, sigmavec);

    /* Calculating sigma --- sigma on rhs is "dT/dp" */
    DMDAVecGetArrayRead (da, Tvec, &T);
    DMDAVecGetArray (da, sigmavec, &sigma);

    LOOP_KJI (
        da,
        sigma[k][j][i] =
            R / p[k] * (R / c_p * T[k][j][i] / p[k] - sigma[k][j][i]);

        if (sigma[k][j][i] < sigmamin) {
        sigma[k][j][i] = sigmamin; });

    DMDAVecRestoreArrayRead (da, Tvec, &T);
    DMDAVecRestoreArray (da, sigmavec, &sigma);

    return (0); }
