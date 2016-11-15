#include "derivatives.h"
#include <strings.h>
#include "petscdmda.h"
#include "context.h"
#include "io.h"

/* *
 * pressure derivative, $ \nabla \times $
 *
 * - one-sided derivatives at the top and bottom boundaries
 */

#undef __FUNCT__
#define __FUNCT__ "pder"
extern PetscErrorCode pder(const PetscScalar hz,const Vec fvec,Vec dvec)
{
        DM                da;
        PetscScalar       wz = 0.5 / hz;
        Vec               avec;
        PetscInt          i,j,k,zs,ys,xs,zm,ym,xm,mz;
        PetscScalar    ***a;
        PetscScalar    ***d;
        PetscErrorCode    ierr;

        PetscFunctionBeginUser;

        ierr = VecGetDM(fvec,&da);CHKERRQ(ierr);
        ierr = DMGetLocalVector(da,&avec);CHKERRQ(ierr);

        ierr = DMGlobalToLocalBegin(da,fvec,INSERT_VALUES,avec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(da,fvec,INSERT_VALUES,avec);CHKERRQ(ierr);

        ierr = DMDAVecGetArrayRead(da,avec,&a);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da,dvec,&d);CHKERRQ(ierr);

        ierr = DMDAGetInfo(da,0,0,0,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
        for (k=zs; k<zs+zm; k++) {
                int k0,k1; PetscScalar wzk;
                if (k==0) {
                        wzk = 2.0 * wz; k1 = k+1; k0 = 0;
                } else if (k==mz-1) {
                        wzk = 2.0 * wz; k1 = k; k0 = k-1;
                } else {
                        wzk = wz; k1 = k+1; k0 = k-1;
                }
                for (j=ys; j<ys+ym; j++) {
                        for (i=xs; i<xs+xm; i++) {
                                d[k][j][i] =
                                        (a[k1][j][i]-a[k0][j][i]) * wzk;
                        }
                }
        }

        ierr = DMDAVecRestoreArrayRead(da,avec,&a);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(da,dvec,&d);CHKERRQ(ierr);

        ierr = DMRestoreLocalVector(da,&avec);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}


/* *
 * Horizontal rotor operator, $ \nabla \times $
 *
 * - one-sided derivatives at the north and south boundaries
 */

#undef __FUNCT__
#define __FUNCT__ "horizontal_rotor"
extern PetscErrorCode horizontal_rotor(Vec Vvec,Vec bvec,Context ctx)
{
        DM             da = ctx->da, da2 = ctx->da2;
        PetscScalar    wx = 0.5 / ctx->hx, wy = 0.5 / ctx->hy;
        PetscInt       my = ctx->my;
        Vec            avec;
        PetscInt       i,j,k,zs,ys,xs,zm,ym,xm;
        PetscScalar    ****a,***b;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;

        ierr = DMGetLocalVector(da2,&avec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalBegin(da2,Vvec,INSERT_VALUES,avec);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(da2,Vvec,INSERT_VALUES,avec);CHKERRQ(ierr);

        ierr = DMDAVecGetArrayDOFRead(da2,avec,&a);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da,bvec,&b);CHKERRQ(ierr);

        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
        for (k=zs; k<zs+zm; k++) {
                for (j=ys; j<ys+ym; j++) {
                        int j0,j1; PetscScalar wyj;
                        if (j==0) {
                                wyj = 2.0 * wy; j1 = j+1; j0 = j;
                        } else if (j==my-1) {
                                wyj = 2.0 * wy; j1 = j; j0 = j-1;
                        } else {
                                wyj = wy; j1 = j+1; j0 = j-1;
                        }
                        for (i=xs; i<xs+xm; i++) {
                                b[k][j][i] =
                                        (a[k][j][i+1][1]-a[k][j][i-1][1]) * wx
                                        - (a[k][j1][i][0]-a[k][j0][i][0]) * wyj;
                        }
                }
        }

        ierr = DMDAVecRestoreArrayDOFRead(da,avec,&a);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(da,bvec,&b);CHKERRQ(ierr);

        ierr = DMRestoreLocalVector(da,&avec);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}


/* *
 * Horizontal advection operator, $ \mathbf{V} \cdot \nabla $
 *
 * - replaces the input field with the result
 * - one-sided derivatives at the north and south boundaries
 */

#undef __FUNCT__
#define __FUNCT__ "horizontal_advection"
        extern PetscErrorCode horizontal_advection(Vec Vvec,Vec bvec,Context ctx)
        {
                DM             da = ctx->da, da2 = ctx->da2;
                PetscScalar    wx = 0.5 / ctx->hx, wy = 0.5 / ctx->hy;
                PetscInt       my = ctx->my;
                Vec            avec;
                PetscInt       i,j,k,zs,ys,xs,zm,ym,xm;
                PetscScalar    ****V,***a,***b;
                PetscErrorCode ierr;

                PetscFunctionBeginUser;

                ierr = DMGetLocalVector(da,&avec);CHKERRQ(ierr);
                ierr = DMGlobalToLocalBegin(da,bvec,INSERT_VALUES,avec);CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(da,bvec,INSERT_VALUES,avec);CHKERRQ(ierr);

                ierr = DMDAVecGetArrayRead(da,avec,&a);CHKERRQ(ierr);
                ierr = DMDAVecGetArrayDOFRead(da2,Vvec,&V);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,bvec,&b);CHKERRQ(ierr);

                ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
                for (k=zs; k<zs+zm; k++) {
                        for (j=ys; j<ys+ym; j++) {
                                int j0,j1; PetscScalar wyj;
                                if (j==0) {
                                        wyj = 2.0 * wy; j1 = j+1; j0 = 0;
                                } else if (j==my-1) {
                                        wyj = 2.0 * wy; j1 = j; j0 = j-1;
                                } else {
                                        wyj = wy; j1 = j+1; j0 = j-1;
                                }
                                for (i=xs; i<xs+xm; i++) {
                                        b[k][j][i]
                                                = V[k][j][i][0] *
                                                (a[k][j][i+1] - a[k][j][i-1]) * wx
                                                + V[k][j][i][1] *
                                                (a[k][j1][i] - a[k][j0][i]) * wyj;
                                }
                        }
                }

                ierr = DMDAVecRestoreArrayRead(da,avec,&a);CHKERRQ(ierr);
                ierr = DMDAVecRestoreArrayDOFRead(da,Vvec,&V);CHKERRQ(ierr);
                ierr = DMDAVecRestoreArray(da,bvec,&b);CHKERRQ(ierr);

                ierr = DMRestoreLocalVector(da,&avec);CHKERRQ(ierr);

                PetscFunctionReturn(0);
        }


/* *
 * Operator coriolis parameter times pressure derivative,
 * $ f \frac{\partial}{\partial p} $
 *
 * - one-sided derivatives at top and bottom boundaries
 */

#undef __FUNCT__
#define __FUNCT__ "fpder"
        extern PetscErrorCode fpder(Vec bvec,Context ctx)
        {
                DM             da = ctx->da;
                PetscInt       mz = ctx->mz;
                PetscScalar    wz = 0.5 / ctx->hz;
                PetscScalar    *f = ctx->Coriolis_parameter;
                Vec            avec;
                PetscInt       i,j,k,zs,ys,xs,zm,ym,xm;
                PetscScalar    ***a,***b;
                PetscErrorCode ierr;

                PetscFunctionBeginUser;

                ierr = DMGetLocalVector(da,&avec);CHKERRQ(ierr);
                ierr = DMGlobalToLocalBegin(da,bvec,INSERT_VALUES,avec);CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(da,bvec,INSERT_VALUES,avec);CHKERRQ(ierr);

                ierr = DMDAVecGetArrayRead(da,avec,&a);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,bvec,&b);CHKERRQ(ierr);

                ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
                for (k=zs; k<zs+zm; k++) {
                        int k0,k1; PetscScalar w;
                        if (k==0) {
                                w = 2.0 * wz; k1 = k+1; k0 = 0;
                        } else if (k==mz-1) {
                                w = 2.0 * wz; k1 = k; k0 = k-1;
                        } else {
                                w = wz; k1 = k+1; k0 = k-1;
                        }
                        for (j=ys; j<ys+ym; j++) {
                                for (i=xs; i<xs+xm; i++) {
                                        b[k][j][i] = f[j]
                                                * (a[k1][j][i] - a[k0][j][i]) * w;
                                }
                        }
                }
                ierr = DMDAVecRestoreArrayRead(da,avec,&a);CHKERRQ(ierr);
                ierr = DMDAVecRestoreArray(da,bvec,&b);CHKERRQ(ierr);

                ierr = DMRestoreLocalVector(da,&avec);CHKERRQ(ierr);

                PetscFunctionReturn(0);
        }
