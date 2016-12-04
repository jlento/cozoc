#include "omegaQG.h"
#include "petscksp.h"
#include "petscdmda.h"
#include <math.h>
#include "context.h"
#include "arrays.h"
#include "constants.h"
#include "io.h"
#include "vecops.h"
#include "derivatives.h"

static PetscErrorCode diff1d(const int n,PetscScalar *x,PetscScalar *f,
                             PetscScalar *dfdx)
{
        PetscFunctionBeginUser;

        dfdx[0] = (f[1] - f[0]) / (x[1] - x[0]);
        for (int i = 1; i < n-1; i++)
                dfdx[i] = (f[i+1] - f[i-1]) / (x[i+1] - x[i-1]);
        dfdx[n-1] = (f[n-1] - f[n-2]) / (x[n-1] - x[n-2]);

        PetscFunctionReturn(0);
}


/*
 * L_{QG}(\omega) = \sigma_0 \nabla^2 \omega
 *                  + f^2 \frac{\partial^2 \omega }{\partial p^2 }
 *
 * \sigma_0 =
 *    \frac{R}{p} \left( \frac{R\overline{T}}{c_pp}
 *                        - \frac{\partial \overline{T}}{\partial p}
 *                \right)
 */


#undef __FUNCT__
#define __FUNCT__ "compute_operator_omega_QG"
extern PetscErrorCode compute_operator_omega_QG(KSP ksp, Mat dummy,
                                                Mat L, void *ctx_p)
{
        Context ctx = (Context)ctx_p;
        PetscScalar    *T_ave,*dTdp,*sigma0,*f2;
        PetscScalar    *p = ctx->Pressure,*f;
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
        PetscScalar    v[7],HyHzdHx,HxHzdHy,HxHydHz;
        const int      sm = 7;
        MatStencil     row,col[sm];
        DM             da;
        //MatNullSpace   nullspace;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;

        ierr = KSPGetDM(ksp, &da); CHKERRQ(ierr);
        ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);
        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);

        ierr = VecGetArray(ctx->f,&f); CHKERRQ(ierr);

        ierr = PetscMalloc4(mz,&T_ave,
                            mz,&dTdp,
                            mz,&sigma0,
                            my,&f2);

        ierr = horizontal_average(ctx,ctx->T,T_ave);

        ierr = diff1d(ctx->mz,p,T_ave,dTdp);

        for (int k=zs; k<zs+zm; k++)
                sigma0[k] = R / p[k] * ( R / c_p * T_ave[k] / p[k] - dTdp[k]);

        for (int j=ys; j<ys+ym; j++)
                f2[j] = f[j] * f[j];

        HyHzdHx = ctx->hy * ctx->hz / ctx->hx;
        HxHzdHy = ctx->hx * ctx->hz / ctx->hy;
        HxHydHz = ctx->hx * ctx->hy / ctx->hz;

        for (k=zs; k<zs+zm; k++) {
                row.k = k;
                col[0].k = k;
                col[1].k = k;
                col[2].k = k;
                col[3].k = k;
                col[4].k = k;
                col[5].k = k-1;
                col[6].k = k+1;
                for (j=ys; j<ys+ym; j++) {
                        row.j = j;
                        col[0].j = j;
                        if(k==0 || k==mz-1 || j==0 || j==my-1) {
                                for (i=xs; i<xs+xm; i++) {
                                        row.i = i;
                                        col[0].i = i;
                                        v[0] = ctx->hx*ctx->hy*ctx->hz;
                                        ierr = MatSetValuesStencil(
                                                L,1,&row,1,col,v,INSERT_VALUES);
                                }
                        } else {
                                col[1].j = j;
                                col[2].j = j;
                                col[3].j = j-1;
                                col[4].j = j+1;
                                col[5].j = j;
                                col[6].j = j;
                                v[0] =  - 2 * sigma0[k] * HyHzdHx
                                        - 2 * sigma0[k] * HxHzdHy
                                        - 2 * f2[j]     * HxHydHz;
                                v[1] = sigma0[k] * HyHzdHx;
                                v[2] = sigma0[k] * HyHzdHx;
                                v[3] = sigma0[k] * HxHzdHy;
                                v[4] = sigma0[k] * HxHzdHy;
                                v[5] = f2[j] * HxHydHz;
                                v[6] = f2[j] * HxHydHz;
                                for (i=xs; i<xs+xm; i++) {
                                        row.i = i;
                                        col[0].i = i;
                                        //col[1].i = (i + mx - 1) % mx ;
                                        //col[2].i = (i + 1) % mx ;
                                        col[1].i = i - 1;
                                        col[2].i = i + 1;
                                        col[3].i = i;
                                        col[4].i = i;
                                        col[5].i = i;
                                        col[6].i = i;
                                        ierr = MatSetValuesStencil(
                                                L,1,&row,7,col,v,INSERT_VALUES);
                                        CHKERRQ(ierr);
                                }
                        }
                }
        }

        ierr = MatAssemblyBegin(L,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(L,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
/* ???
   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
   ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
   ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
*/
        ierr = VecRestoreArray(ctx->f,&f); CHKERRQ(ierr);
        ierr = PetscFree4(T_ave,dTdp,sigma0,f2); CHKERRQ(ierr);

        PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "geostrophic_wind_v"
extern PetscErrorCode geostrophic_wind_v(Vec Z,Vec f,PetscScalar hx,Vec v)
{
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;
        ierr = xder(Z,hx,v);CHKERRQ(ierr);
        ierr = ydivide(v,v,f);CHKERRQ(ierr);
        ierr = VecScale(v,g);CHKERRQ(ierr);
        PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "geostrophic_wind_u"
        extern PetscErrorCode geostrophic_wind_u(Context ctx,Vec u)
        {
                DM da;
                PetscInt i,j,k,zs,ys,xs,zm,ym,xm;
                PetscScalar *f,w;
                PetscScalar ***ua;
                PetscErrorCode  ierr;

                PetscFunctionBeginUser;

                ierr = VecGetDM(u,&da);
                ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

                ierr = yder(ctx->Z,ctx->hy,u);CHKERRQ(ierr);

                ierr = VecGetArray(ctx->f,&f);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,u,&ua);CHKERRQ(ierr);
                for (k=zs; k<zs+zm; k++) {
                        for (j=ys; j<ys+ym; j++) {
                                w = -g / f[j];
                                for (i=xs; i<xs+xm; i++) {
                                        ua[k][j][i] *= w;
                                }
                        }
                }
                ierr = DMDAVecRestoreArray(da,u,&ua);CHKERRQ(ierr);
                ierr = VecRestoreArray(ctx->f,&f);CHKERRQ(ierr);
                PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "surface_factor"
        extern PetscErrorCode surface_factor(Context ctx,Vec s)
        {
                PetscInt i,j,k,zs,ys,xs,zm,ym,xm;
                PetscScalar ***sa,*p=ctx->Pressure,**psfc;
                PetscErrorCode  ierr;

                PetscFunctionBeginUser;

                ierr = DMDAGetCorners(ctx->da,&xs,&ys,&zs,&xm,&ym,&zm);
                CHKERRQ(ierr);

                ierr = DMDAVecGetArray(ctx->da,s,&sa);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(ctx->daxy,ctx->psfc,&psfc);CHKERRQ(ierr);
                for (k=1; k<ctx->mz-1; k++) {
                        for (j=ys; j<ys+ym; j++) {
                                for (i=xs; i<xs+xm; i++) {
                                        if (psfc[j][i] <= (p[k]+p[k+1])/2 ) {
                                                sa[k][j][i] = 0.0;
                                        } else if (psfc[j][i] <= (p[k]+p[k-1])/2) {
                                                sa[k][j][i] =
                                                        (p[k]+p[k+1]-2.0*psfc[j][i])
                                                        / (p[k+1]-p[k-1]);
                                        } else {
                                                sa[k][j][i] = 1.0;
                                        }
                                }
                        }
                }
                k = ctx->mz-1;
                for (j=ys; j<ys+ym; j++) {
                        for (i=xs; i<xs+xm; i++)
                                sa[k][j][i] = 1.0;
                }
                k = 0;
                for (j=ys; j<ys+ym; j++) {
                        for (i=xs; i<xs+xm; i++) {
                                if (psfc[j][i] > (p[k]+p[k+1])/2) {
                                        sa[k][j][i] =
                                                (p[k]+p[k+1]-2.0*psfc[j][i])
                                                / (p[k+1]-p[k]) / 2.0;
                                } else {
                                        sa[k][j][i] = 1.0;
                                }
                        }
                }
                ierr = DMDAVecRestoreArray(ctx->daxy,ctx->psfc,&psfc);
                CHKERRQ(ierr);
                ierr = DMDAVecRestoreArray(ctx->da,s,&sa);CHKERRQ(ierr);
                PetscFunctionReturn(0);
        }


#undef __FUNCT__
#define __FUNCT__ "compute_vorticity_advection"
extern PetscErrorCode compute_vorticity_advection(Context ctx,Vec u,Vec v,
                                                  Vec F_V)
{
        DM da;
        Vec s;
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;

        ierr = VecGetDM(F_V,&da);CHKERRQ(ierr);

        ierr = VecDuplicate(F_V,&s);CHKERRQ(ierr);

        /* F_V = zeta_g */
        ierr = curl(u,v,ctx->hx,ctx->hy,F_V);CHKERRQ(ierr);
        /* F_V = zeta_g + f */
        ierr = Vec2Vec123AXPY(F_V,1.0,ctx->f);CHKERRQ(ierr);
        /* s = adv := V_G . grad(zeta_g + f) */
        ierr = advection(u,v,F_V,ctx->hx,ctx->hy,F_V);CHKERRQ(ierr);
        ierr = surface_factor(ctx,s);CHKERRQ(ierr);
        ierr = VecPointwiseMult(F_V,s,F_V);CHKERRQ(ierr);
        /* F_V = f * d/dp adv */
        ierr = fpder(F_V,ctx);CHKERRQ(ierr);

        ierr = VecDestroy(&s);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "compute_rhs_omega_QG"
        extern PetscErrorCode compute_rhs_omega_QG(KSP ksp,Vec b,void *ctx_p)
        {
                Context ctx = (Context)ctx_p;
                Vec Z = ctx->Z, f = ctx->f;
                PetscScalar hx = ctx->hx;

                Vec F_V, F_T;

                Vec u,v;
                PetscErrorCode  ierr;

                PetscFunctionBeginUser;

                F_V = b;
                ierr = VecDuplicate(b,&F_T);CHKERRQ(ierr);
                ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
                ierr = VecDuplicate(b,&v);CHKERRQ(ierr);

                ierr = geostrophic_wind_u(ctx,u);CHKERRQ(ierr);
                ierr = geostrophic_wind_v(Z,f,hx,v);CHKERRQ(ierr);

                ierr = compute_vorticity_advection(ctx,u,v,F_V);CHKERRQ(ierr);
                ierr = VecScale(F_V,ctx->hx*ctx->hy*ctx->hz);CHKERRQ(ierr);

//ierr = compute_thermal_advection(ctx,u,v,F_T);CHKERRQ(ierr);

//ierr = VecAXPY(F_V,1.0,F_T);CHKERRQ(ierr);

                ierr = VecDestroy(&v);CHKERRQ(ierr);
                ierr = VecDestroy(&u);CHKERRQ(ierr);
                ierr = VecDestroy(&F_T);CHKERRQ(ierr);
                PetscFunctionReturn(0);
        }
