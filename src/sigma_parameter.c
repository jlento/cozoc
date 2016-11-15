#include "sigma_parameter.h"
#include "petscdmda.h"
#include "derivatives.h"
#include "loops.h"
#include "constants.h"

#undef __FUNCT__
#define __FUNCT__ "sigma_parameter"
extern PetscErrorCode sigma_parameter(Context ctx)
{
        DM                da       = ctx->da;
        PetscScalar       hz       = ctx->hz;
        PetscScalar      *p        = ctx->Pressure;
        Vec               Tvec     = ctx->T;
        Vec               sigmavec = ctx->Sigma_parameter;
        PetscScalar    ***T;
        PetscScalar    ***sigma;
        PetscErrorCode    ierr;

        PetscFunctionBeginUser;

        /* first sigma holds auxiliary variable dT/dp, temporarily */
        ierr = pder(hz,Tvec,sigmavec);CHKERRQ(ierr);

        /* Calculating sigma --- sigma on rhs is "dT/dp" */
        ierr = DMDAVecGetArrayRead(da,Tvec,&T);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da,sigmavec,&sigma);CHKERRQ(ierr);
        LOOP_KJI(da,
                 sigma[k][j][i] =
                 R / p[k] * ( R / c_p * T[k][j][i] / p[k] - sigma[k][j][i]);
                 if (sigma[k][j][i] < sigmamin) {
                         sigma[k][j][i] = sigmamin;
                 }
                );
        ierr = DMDAVecRestoreArrayRead(da,Tvec,&T);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(da,sigmavec,&sigma);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}
