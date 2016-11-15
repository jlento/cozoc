#include "vorticity.h"
#include "petscdmda.h"
#include "derivatives.h"
#include "loops.h"
#include "constants.h"

#undef __FUNCT__
#define __FUNCT__ "vorticity"
extern PetscErrorCode vorticity(Context ctx)
{
        DM                da       = ctx->da;
        DM                da2      = ctx->da2;
        Vec               Vvec     = ctx->Horizontal_wind;
        Vec               sigmavec = ctx->Sigma_parameter;
        Vec               zetavec  = ctx->Vorticity;
        PetscScalar       hz       = ctx->hz;
        size_t            mz       = ctx->mz;
        PetscScalar      *f        = ctx->Coriolis_parameter;
        PetscScalar    ***sigma;
        PetscScalar   ****V;
        PetscScalar    ***zeta;
        PetscScalar       zetap,dudp,dvdp,dVdp2;
        PetscErrorCode    ierr;

        PetscFunctionBeginUser;

        ierr = horizontal_rotor(Vvec,zetavec,ctx);CHKERRQ(ierr);

        ierr = DMDAVecGetArray(da,zetavec,&zeta);CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(da,sigmavec,&sigma);CHKERRQ(ierr);
        ierr = DMDAVecGetArrayDOFRead(da2,Vvec,&V);CHKERRQ(ierr);
        LOOP_KJI(da,
                 if (f[j] > 1.0e-7) {
                         int k0; int k1; PetscScalar w;
                         if (k == 0) {
                                 k1 = 1; k0 = 0; w = 1.0 / (hz*hz);
                         } else if (k == mz-1) {
                                 k1 = mz-1; k0 = mz-2; w = 1.0 / (hz*hz);
                         } else {
                                 k1 = k+1; k0 = k-1; w = 0.25 / (hz*hz);
                         }
                         dudp = w*(V[k1][j][i][0] - V[k0][j][i][0]);
                         dvdp = w*(V[k1][j][i][1] - V[k0][j][i][1]);
                         dVdp2 = dudp*dudp + dvdp*dvdp;
                         zetap = etamin
                                 + f[j] / (4.0 * sigma[k][j][i]) * dVdp2
                                 - f[j];
                         if (zetap > zeta[k][j][i])
                                 zeta[k][j][i] = zetap;
                 }
                );
        ierr = DMDAVecRestoreArray(da,zetavec,&zeta);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(da,sigmavec,&sigma);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayDOFRead(da2,Vvec,&V);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}
