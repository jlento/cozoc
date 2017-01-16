#include "constants.h"
#include "context.h"
#include "ops.h"
#include "omegaQG.h"
#include "petscdmda.h"
#include "petscksp.h"

/*
 * L_{QG}(\omega) = \sigma_0 \nabla^2 \omega
 *                  + f^2 \frac{\partial^2 \omega }{\partial p^2 }
 *
 * \sigma_0 =
 *    \frac{R}{p} \left( \frac{R\overline{T}}{c_pp}
 *                        - \frac{\partial \overline{T}}{\partial p}
 *                \right)
 */

extern PetscErrorCode omega_qg_compute_operator (
    KSP ksp, Mat dummy, Mat L, void* ctx_p) {

    Context      ctx = (Context) ctx_p;
    PetscScalar* p   = ctx->Pressure;
    PetscScalar* f   = ctx->Coriolis_parameter;
    Vec          T   = ctx->Temperature;
    const double R   = Specific_gas_constant_of_dry_air;
    const double c_p = Specific_heat_of_dry_air;
    PetscScalar *T_ave, *dTdp, *sigma0, *f2;
    PetscInt     i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
    PetscScalar  v[7], HyHzdHx, HxHzdHy, HxHydHz;
    const int    sm = 7;
    MatStencil   row, col[sm];
    DM           da;
    // MatNullSpace   nullspace;

    KSPGetDM (ksp, &da);
    DMDAGetInfo (da, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    PetscMalloc4 (mz, &T_ave, mz, &dTdp, mz, &sigma0, my, &f2);

    horizontal_average (ctx, T, T_ave);

    diff1d (ctx->mz, p, T_ave, dTdp);

    for (int k    = zs; k < zs + zm; k++)
        sigma0[k] = R / p[k] * (R / c_p * T_ave[k] / p[k] - dTdp[k]);

    for (int j = ys; j < ys + ym; j++)
        f2[j]  = f[j] * f[j];

    HyHzdHx = ctx->hy * ctx->hz / ctx->hx;
    HxHzdHy = ctx->hx * ctx->hz / ctx->hy;
    HxHydHz = ctx->hx * ctx->hy / ctx->hz;

    for (k = zs; k < zs + zm; k++) {
        row.k    = k;
        col[0].k = k;
        col[1].k = k;
        col[2].k = k;
        col[3].k = k;
        col[4].k = k;
        col[5].k = k - 1;
        col[6].k = k + 1;

        for (j = ys; j < ys + ym; j++) {
            row.j    = j;
            col[0].j = j;

            if (k == 0 || k == mz - 1 || j == 0 || j == my - 1) {
                for (i = xs; i < xs + xm; i++) {
                    row.i    = i;
                    col[0].i = i;
                    v[0]     = ctx->hx * ctx->hy * ctx->hz;
                    MatSetValuesStencil (
                        L, 1, &row, 1, col, v, INSERT_VALUES); } }
            else {
                col[1].j = j;
                col[2].j = j;
                col[3].j = j - 1;
                col[4].j = j + 1;
                col[5].j = j;
                col[6].j = j;
                v[0]     = -2 * sigma0[k] * HyHzdHx -
                           2 * sigma0[k] * HxHzdHy - 2 * f2[j] * HxHydHz;
                v[1] = sigma0[k] * HyHzdHx;
                v[2] = sigma0[k] * HyHzdHx;
                v[3] = sigma0[k] * HxHzdHy;
                v[4] = sigma0[k] * HxHzdHy;
                v[5] = f2[j] * HxHydHz;
                v[6] = f2[j] * HxHydHz;

                for (i = xs; i < xs + xm; i++) {
                    row.i    = i;
                    col[0].i = i;
                    // col[1].i = (i + mx - 1) % mx ;
                    // col[2].i = (i + 1) % mx ;
                    col[1].i = i - 1;
                    col[2].i = i + 1;
                    col[3].i = i;
                    col[4].i = i;
                    col[5].i = i;
                    col[6].i = i;
                    MatSetValuesStencil (
                        L, 1, &row, 7, col, v, INSERT_VALUES); } } } }

    MatAssemblyBegin (L, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd (L, MAT_FINAL_ASSEMBLY);
    /* ???
       ierr =
       MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
       MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
       MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
    */

    PetscFree4 (T_ave, dTdp, sigma0, f2);

    return (0); }


/* *
 * Geostrophic wind
 *
 * V_G = \left(
 *              -\frac{g}{f}\frac{\partial Z}{\partial y},
 *               \frac{g}{f}\frac{\partial Z}{\partial x}
 *       \right)
 */

extern int geostrophic_wind (Vec Vvec, Context ctx) {

    DM           da   = ctx->da;
    DM           da2  = ctx->da2;
    Vec          Zvec = ctx->Geopotential_height;
    PetscScalar* f    = ctx->Coriolis_parameter;
    PetscScalar  wx   = 0.5 / ctx->hx;
    PetscScalar  wy   = 0.5 / ctx->hy;
    PetscInt     my   = ctx->my;
    const double g    = Gravitational_acceleration;

    Vec            Zloc;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ***Z, ****V;

    DMGetLocalVector (da, &Zloc);
    DMGlobalToLocalBegin (da, Zvec, INSERT_VALUES, Zloc);
    DMGlobalToLocalEnd (da, Zvec, INSERT_VALUES, Zloc);

    DMDAVecGetArrayRead (da, Zloc, &Z);
    DMDAVecGetArrayDOF (da2, Vvec, &V);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    int         j0, j1;
    PetscScalar wyj, wxi;

    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            if (j == 0) {
                wyj = -g / f[j] * 2.0 * wy;
                j1  = j + 1;
                j0  = j; }
            else if (j == my - 1) {
                wyj = -g / f[j] * 2.0 * wy;
                j1  = j;
                j0  = j - 1; }
            else {
                wyj = -g / f[j] * wy;
                j1  = j + 1;
                j0  = j - 1; }

            wxi = g / f[j] * wx;

            for (i = xs; i < xs + xm; i++) {
                V[k][j][i][0] = wyj * (Z[k][j1][i] - Z[k][j0][i]);
                V[k][j][i][1] = wxi * (Z[k][j][i + 1] - Z[k][j][i - 1]); } } }


    DMDAVecRestoreArrayRead (da, Zloc, &Z);
    DMDAVecRestoreArrayDOF (da, Vvec, &V);

    DMRestoreLocalVector (da, &Zloc);

    return (0); }


/* Geostrophic vorticity
 *
 * zeta_g = \nabla \rot V_g
 *        = \left(
 *                   \frac{\partial v_g}{\partial x}
 *                 - \frac{\partial u_g}{\partial y}
            \right)
 */

extern PetscErrorCode geostrophic_vorticity (
    Vec result, Vec V_g, Context ctx) {
    DM           da  = ctx->da;
    DM           da2 = ctx->da2;
    size_t       my  = ctx->my;
    PetscScalar  hx  = ctx->hx;
    PetscScalar  hy  = ctx->hy;

    horizontal_rotor (da, da2, my, hx, hy, V_g, result);

    return (0); }


extern PetscErrorCode surface_factor (Context ctx, Vec s) {

    DM           da      = ctx->da;
    DM           daxy    = ctx->daxy;
    PetscScalar* p       = ctx->Pressure;
    Vec          PSFCVec = ctx->Surface_pressure;

    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ***sa, **psfc;

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    DMDAVecGetArray (da, s, &sa);
    DMDAVecGetArray (daxy, PSFCVec, &psfc);

    for (k = 1; k < (int) ctx->mz - 1; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                if (psfc[j][i] <= (p[k] + p[k + 1]) / 2) {
                    sa[k][j][i] = 0.0; }
                else if (psfc[j][i] <= (p[k] + p[k - 1]) / 2) {
                    sa[k][j][i] = (p[k] + p[k + 1] - 2.0 * psfc[j][i]) /
                                  (p[k + 1] - p[k - 1]); }
                else {
                    sa[k][j][i] = 1.0; } } } }

    k = ctx->mz - 1;

    for (j = ys; j < ys + ym; j++) {
        for (i          = xs; i < xs + xm; i++)
            sa[k][j][i] = 1.0; }

    k = 0;

    for (j = ys; j < ys + ym; j++) {
        for (i = xs; i < xs + xm; i++) {
            if (psfc[j][i] > (p[k] + p[k + 1]) / 2) {
                sa[k][j][i] = (p[k] + p[k + 1] - 2.0 * psfc[j][i]) /
                              (p[k + 1] - p[k]) / 2.0; }
            else {
                sa[k][j][i] = 1.0; } } }

    DMDAVecRestoreArray (daxy, PSFCVec, &psfc);
    DMDAVecRestoreArray (da, s, &sa);
    return (0); }


/* *
 * F_{V(QG)} = f \frac{\partial}{\partial p}
 *               \left[
 *                       \mathbf{V_g} \cdot \nabla \left( \zeta_g +f
 * \right)
 *               \right]
 *
 *
 * F_{T(QG)} = \frac{R}{p} \nabla^2 ( V_g \cdot \nabla T )
 *
 */

extern PetscErrorCode omega_qg_compute_rhs (
    KSP ksp, Vec F, void* ctx_p) {

    Context      ctx = (Context) ctx_p;
    DM           da2 = ctx->da2;
    DM           da  = ctx->da;
    size_t       mz      = ctx->mz;
    PetscScalar* p       = ctx->Pressure;
    Vec          T   = ctx->Temperature;
    PetscScalar* f   = ctx->Coriolis_parameter;
    Vec          V_g;
    Vec          F_T;

    DMGetGlobalVector (da2, &V_g);
    DMGetGlobalVector (da, &F_T);

    geostrophic_wind (V_g, ctx);

    /* Vorticity advection forcing F_V*/
    geostrophic_vorticity (F, V_g, ctx);    // F = zeta_g
    field_array1d_add (F, f, DMDA_Y);       // (f+)
    horizontal_advection (F, V_g, ctx);     // (V_g \cdot \nabla)
    fpder (da, mz, f, p, F);                // (f * d/dp)

    /* Temperature advection forcing F_T*/
    VecCopy (T, F_T);                        // F_T = T
    horizontal_advection (F_T, V_g, ctx);    // (V_g \cdot \nabla)
    plaplace (F_T, ctx);                     // (R/p * \nabla^2)

    /* F = F_V + F_T */
    VecAXPY (F, 1.0, F_T);
    VecScale (F, ctx->hx * ctx->hy * ctx->hz);

    DMRestoreGlobalVector (ctx->da2, &V_g);
    DMRestoreGlobalVector (ctx->da, &F_T);
    return (0); }
