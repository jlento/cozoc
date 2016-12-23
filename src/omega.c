#include "context.h"
//#include "coriolis_parameter.h"
#include "derivatives.h"
#include "field.h"
#include "io.h"
#include "omega.h"
#include <petscdmda.h>
#include <petscksp.h>

#include "omega_stencil.inc"


char* omega_component_id_string[N_OMEGA_COMPONENTS] = {"ome_v",
                                                       "ome_t" };

PetscErrorCode (*omega_compute_rhs[N_OMEGA_COMPONENTS]) (
    KSP ksp, Vec b, void* ctx_p) = {omega_compute_rhs_F_V,
                                    omega_compute_rhs_F_T };


/* *
 * The LHS of the generalized omega equation
 *
 */

extern PetscErrorCode omega_compute_operator (
    KSP ksp, Mat A, Mat B, void* ctx_p) {

    Context              ctx = (Context) ctx_p;
    DM                   da  = ctx->da;
    DM                   da2 = ctx->da2;
    PetscScalar          hx  = ctx->hx;
    PetscScalar          hy  = ctx->hy;
    PetscScalar          hz  = ctx->hz;
    PetscScalar*         f   = ctx->Coriolis_parameter;
    Vec                  sigmavec, zetavec, Vvec;
    const PetscScalar ***sigma, ***zeta, ****V;
    PetscInt             my = ctx->my;
    PetscInt             mz = ctx->mz;
    PetscInt             n;
    MatStencil           row, col[15];
    PetscScalar          w[15];
    PetscInt             xs, ys, zs, xm, ym, zm;

    DMGetLocalVector (da, &sigmavec);
    DMGetLocalVector (da, &zetavec);
    DMGetLocalVector (da2, &Vvec);

    DMGlobalToLocalBegin (
        da, ctx->Sigma_parameter, INSERT_VALUES, sigmavec);
    DMGlobalToLocalEnd (
        da, ctx->Sigma_parameter, INSERT_VALUES, sigmavec);
    DMGlobalToLocalBegin (da, ctx->Vorticity, INSERT_VALUES, zetavec);
    DMGlobalToLocalEnd (da, ctx->Vorticity, INSERT_VALUES, zetavec);
    DMGlobalToLocalBegin (
        da2, ctx->Horizontal_wind, INSERT_VALUES, Vvec);
    DMGlobalToLocalEnd (da2, ctx->Horizontal_wind, INSERT_VALUES, Vvec);

    DMDAVecGetArrayRead (da, sigmavec, &sigma);
    DMDAVecGetArrayRead (da, zetavec, &zeta);
    DMDAVecGetArrayDOFRead (da2, Vvec, &V);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                row.i = i;
                row.j = j;
                row.k = k;

                if (k == 0 || k == mz - 1 || j == 0 || j == my - 1) {
                    w[0] = hx * hy * hz;
                    MatSetValuesStencil (
                        B, 1, &row, 1, &row, w, INSERT_VALUES); }
                else {
                    omega_stencil (
                        i,
                        j,
                        k,
                        hx,
                        hy,
                        hz,
                        f,
                        sigma,
                        zeta,
                        V,
                        &n,
                        col,
                        w);
                    MatSetValuesStencil (
                        B, 1, &row, n, col, w, INSERT_VALUES); } } } }


    MatAssemblyBegin (B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd (B, MAT_FINAL_ASSEMBLY);

    DMDAVecRestoreArrayRead (da, sigmavec, &sigma);
    DMDAVecRestoreArrayRead (da, zetavec, &zeta);
    DMDAVecRestoreArrayDOFRead (da2, Vvec, &V);

    DMRestoreLocalVector (da, &sigmavec);
    DMRestoreLocalVector (da, &zetavec);
    DMRestoreLocalVector (da2, &Vvec);

    return (0); }


/* *
 * Vorticity advection forcing
 *
 * F_V = f \frac{\partial}{\partial p}
 *               \left[
 *                       \mathbf{V} \cdot \nabla \left( \zeta +f \right)
 *               \right]
 */

extern PetscErrorCode omega_compute_rhs_F_V (
    KSP ksp, Vec b, void* ctx_p) {

    Context      ctx  = (Context) ctx_p;
    Vec          zeta = ctx->Vorticity;
    Vec          V    = ctx->Horizontal_wind;
    PetscScalar* f    = ctx->Coriolis_parameter;
    PetscScalar  hx   = ctx->hx;
    PetscScalar  hy   = ctx->hy;
    PetscScalar  hz   = ctx->hz;

    VecCopy (zeta, b);
    field_array1d_add (b, f, DMDA_Y);
    horizontal_advection (b, V, ctx);
    /* TODO: surface factor */
    fpder (b, ctx);
    VecScale (b, hx * hy * hz);
    /*        write3Ddump(ctx,"b",b); */

    return (0); }


/* *
 * Temperature advection forcing
 *
 * F_T = \frac{R}{p} \nabla^2 ( \mathbf{V} \cdot \nabla \mathbf{T} )
 *
 */

extern PetscErrorCode omega_compute_rhs_F_T (
    KSP ksp, Vec b, void* ctx_p) {

    Context      ctx = (Context) ctx_p;
    Vec          T   = ctx->Temperature;
    Vec          V   = ctx->Horizontal_wind;
    PetscScalar  hx  = ctx->hx;
    PetscScalar  hy  = ctx->hy;
    PetscScalar  hz  = ctx->hz;

    VecCopy (T, b);
    horizontal_advection (b, V, ctx);
    plaplace (b, ctx);
    VecScale (b, hx * hy * hz);

    return (0); }
