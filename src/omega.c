#include "omega.h"
#include "petscksp.h"
#include "petscdmda.h"
#include "context.h"
#include "coriolis_parameter.h"
#include "derivatives.h"
#include "loops.h"
#include "io.h"

#include "omega_stencil.inc"

/* *
 * F_V = f \frac{\partial}{\partial p}
 *               \left[
 *                       \mathbf{V} \cdot \nabla \left( \zeta +f \right)
 *               \right]
 */

#undef __FUNCT__
#define __FUNCT__ "omega_compute_rhs_F_V"
extern PetscErrorCode omega_compute_rhs_F_V(KSP ksp,Vec b,void *ctx_p)
{
        Context        ctx  = (Context)ctx_p;
        Vec            zeta = ctx->Vorticity;
        Vec            V    = ctx->Horizontal_wind;
        PetscScalar    hx   = ctx->hx;
        PetscScalar    hy   = ctx->hy;
        PetscScalar    hz   = ctx->hz;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        ierr = VecCopy(zeta,b);CHKERRQ(ierr);
        ierr = coriolis_parameter_add(b,ctx);CHKERRQ(ierr);
        ierr = horizontal_advection(V,b,ctx);CHKERRQ(ierr);
        /* TODO: surface factor */
        ierr = fpder(b,ctx);CHKERRQ(ierr);
        ierr = VecScale(b,hx*hy*hz);CHKERRQ(ierr);
/*        write3Ddump(ctx,"b",b); */

        PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "omega_compute_operator"
extern PetscErrorCode omega_compute_operator(KSP ksp,Mat A,Mat B,void *ctx_p)
{
        Context              ctx = (Context)ctx_p;
        DM                   da  = ctx->da;
        DM                   da2 = ctx->da2;
        PetscScalar          hx  = ctx->hx;
        PetscScalar          hy  = ctx->hy;
        PetscScalar          hz  = ctx->hz;
        PetscScalar         *f   = ctx->Coriolis_parameter;
        Vec                  sigmavec,zetavec,Vvec;
        const PetscScalar ***sigma,***zeta,****V;
        PetscInt             my = ctx->my;
        PetscInt             mz = ctx->mz;
        PetscInt             n;
        MatStencil           row,col[15];
        PetscScalar          w[15];
        PetscErrorCode       ierr;

        PetscFunctionBeginUser;

        ierr = DMGetLocalVector(da,&sigmavec);CHKERRQ(ierr);
        ierr = DMGetLocalVector(da,&zetavec);CHKERRQ(ierr);
        ierr = DMGetLocalVector(da2,&Vvec);CHKERRQ(ierr);

        ierr = DMGlobalToLocalBegin(da,ctx->Sigma_parameter,INSERT_VALUES,
                                    sigmavec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(da,ctx->Sigma_parameter,INSERT_VALUES,
                                  sigmavec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalBegin(da,ctx->Vorticity,INSERT_VALUES,
                                    zetavec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(da,ctx->Vorticity,INSERT_VALUES,
                                  zetavec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalBegin(da2,ctx->Horizontal_wind,INSERT_VALUES,
                                    Vvec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(da2,ctx->Horizontal_wind,INSERT_VALUES,
                                  Vvec);CHKERRQ(ierr);

        ierr = DMDAVecGetArrayRead(da,sigmavec,&sigma);CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(da,zetavec,&zeta);CHKERRQ(ierr);
        ierr = DMDAVecGetArrayDOFRead(da2,Vvec,&V);CHKERRQ(ierr);

        LOOP_KJI(da,
                 row.i = i; row.j = j; row.k = k;
                 if (k==0 || k==mz-1 || j==0 || j==my-1) {
                         w[0] = hx*hy*hz;
                         ierr = MatSetValuesStencil(
                                 B,1,&row,1,&row,w,
                                 INSERT_VALUES);
                 } else {
                         omega_stencil(
                                 i,j,k,hx,hy,hz,
                                 f,sigma,zeta,V,
                                 &n,col,w);
                         ierr = MatSetValuesStencil(
                                 B,1,&row,n,col,w,
                                 INSERT_VALUES);
                 }
                );

        ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

        ierr = DMDAVecRestoreArrayRead(da,sigmavec,&sigma);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(da,zetavec,&zeta);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayDOFRead(da2,Vvec,&V);CHKERRQ(ierr);

        ierr = DMRestoreLocalVector(da,&sigmavec);CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(da,&zetavec);CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(da2,&Vvec);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}
