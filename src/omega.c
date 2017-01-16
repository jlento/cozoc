#include "constants.h"
#include "context.h"
#include "ops.h"
#include "omega.h"
#include <petscdmda.h>
#include <petscksp.h>

#include "omega_stencil.inc"


char* omega_component_id_string[N_OMEGA_COMPONENTS] = {
    "ome_v", "ome_t", "ome_f", "ome_q", "ome_a" };

PetscErrorCode (*omega_compute_rhs[N_OMEGA_COMPONENTS]) (
    KSP ksp, Vec b, void* ctx_p) = {omega_compute_rhs_F_V,
                                    omega_compute_rhs_F_T,
                                    omega_compute_rhs_F_F,
                                    omega_compute_rhs_F_Q,
                                    omega_compute_rhs_F_A };


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
    DM           da   = ctx->da;
    size_t       mz   = ctx->mz;
    PetscScalar* p    = ctx->Pressure;
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
    fpder (da, mz, f, p, b);
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

    Context     ctx = (Context) ctx_p;
    Vec         T   = ctx->Temperature;
    Vec         V   = ctx->Horizontal_wind;
    PetscScalar hx  = ctx->hx;
    PetscScalar hy  = ctx->hy;
    PetscScalar hz  = ctx->hz;

    VecCopy (T, b);
    horizontal_advection (b, V, ctx);
    plaplace (b, ctx);
    VecScale (b, hx * hy * hz);

    return (0); }


/* *
 * Friction forcing
 *
 * F_F = - f \frac{\partial}{\partial p}
 *               \left(
 *                       \mathbf{k} \cdot \nabla \times \mathbf{F}
 *               \right)
 */

extern PetscErrorCode omega_compute_rhs_F_F (
    KSP ksp, Vec b, void* ctx_p) {

    Context      ctx = (Context) ctx_p;
    DM           da  = ctx->da;
    DM           da2 = ctx->da2;
    size_t       my  = ctx->my;
    size_t       mz  = ctx->mz;
    PetscScalar* p   = ctx->Pressure;
    PetscScalar* f   = ctx->Coriolis_parameter;
    Vec          F   = ctx->Friction;
    PetscScalar  hx  = ctx->hx;
    PetscScalar  hy  = ctx->hy;
    PetscScalar  hz  = ctx->hz;

    horizontal_rotor (da, da2, my, hx, hy, F, b);
    fpder (da, mz, f, p, b);
    VecScale (b, -hx * hy * hz);

    return (0); }


/* *
 * Diabatic heating forcing
 *
 * F_Q = -\frac{R}{c_p p} \nabla^2 ( \mathbf{V} \cdot \nabla \mathbf{T}
 * )
 *
 */

extern PetscErrorCode omega_compute_rhs_F_Q (
    KSP ksp, Vec b, void* ctx_p) {

    Context     ctx = (Context) ctx_p;
    Vec         Q   = ctx->Diabatic_heating;
    PetscScalar hx  = ctx->hx;
    PetscScalar hy  = ctx->hy;
    PetscScalar hz  = ctx->hz;
    PetscScalar c_p = Specific_heat_of_dry_air;

    VecCopy (Q, b);
    plaplace (b, ctx);
    VecScale (b, -hx * hy * hz / c_p);

    return (0); }


/* *
 * Ageostrophic vorticity tendency forcing
 *
 * F_A = \left[
 *           f \frac{\partial}{\partial p} \left(
 *               \frac{\partial \zeta}{\partial t} \right)
 *           + \frac{R}{p} \nabla^2 \frac{\partial T}{\partial t}
 *       \right]
 */

extern PetscErrorCode omega_compute_rhs_F_A (
    KSP ksp, Vec b, void* ctx_p) {

    Context      ctx     = (Context) ctx_p;
    DM           da      = ctx->da;
    size_t       mz      = ctx->mz;
    PetscScalar* p       = ctx->Pressure;
    PetscScalar* f       = ctx->Coriolis_parameter;
    Vec          dzetadt = ctx->Vorticity_tendency;
    Vec          dTdt    = ctx->Temperature_tendency;
    PetscScalar  hx      = ctx->hx;
    PetscScalar  hy      = ctx->hy;
    PetscScalar  hz      = ctx->hz;
    Vec          tmpvec;

    DMGetGlobalVector (da, &tmpvec);
    VecCopy (dTdt, b);
    plaplace (b, ctx);
    VecCopy (dzetadt, tmpvec);

    fpder (da, mz, f, p, tmpvec);
    VecAXPY (b, 1.0, tmpvec);
    DMRestoreGlobalVector (da, &tmpvec);
    VecScale (b, hx * hy * hz);

    return (0); }
