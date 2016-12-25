#include "constants.h"
#include "context.h"
#include "derivatives.h"
#include "io.h"
#include "petscdmda.h"
#include <strings.h>


extern PetscErrorCode diff1d (
    const int n, PetscScalar* x, PetscScalar* f, PetscScalar* dfdx) {

    dfdx[0] = (f[1] - f[0]) / (x[1] - x[0]);

    for (int i  = 1; i < n - 1; i++)
        dfdx[i] = (f[i + 1] - f[i - 1]) / (x[i + 1] - x[i - 1]);

    dfdx[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);

    return (0); }


/* *
 * pressure derivative, $ \nabla \times $
 *
 * - one-sided derivatives at the top and bottom boundaries
 */

extern PetscErrorCode pder (
    const PetscScalar hz, const Vec fvec, Vec dvec) {

    DM             da;
    PetscScalar    wz = 0.5 / hz;
    Vec            avec;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm, mz;
    PetscScalar*** a;
    PetscScalar*** d;

    VecGetDM (fvec, &da);
    DMGetLocalVector (da, &avec);

    DMGlobalToLocalBegin (da, fvec, INSERT_VALUES, avec);
    DMGlobalToLocalEnd (da, fvec, INSERT_VALUES, avec);

    DMDAVecGetArrayRead (da, avec, &a);
    DMDAVecGetArray (da, dvec, &d);

    DMDAGetInfo (da, 0, 0, 0, &mz, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (k = zs; k < zs + zm; k++) {
        int         k0, k1;
        PetscScalar wzk;

        if (k == 0) {
            wzk = 2.0 * wz;
            k1  = k + 1;
            k0  = 0; }
        else if (k == mz - 1) {
            wzk = 2.0 * wz;
            k1  = k;
            k0  = k - 1; }
        else {
            wzk = wz;
            k1  = k + 1;
            k0  = k - 1; }

        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                d[k][j][i] = (a[k1][j][i] - a[k0][j][i]) * wzk; } } }

    DMDAVecRestoreArrayRead (da, avec, &a);
    DMDAVecRestoreArray (da, dvec, &d);

    DMRestoreLocalVector (da, &avec);

    return (0); }


/* *
 * Horizontal rotor operator, $ \nabla \times $
 *
 * - one-sided derivatives at the north and south boundaries
 */

extern PetscErrorCode horizontal_rotor (
    Vec Vvec, Vec bvec, Context ctx) {

    DM              da  = ctx->da;
    DM              da2 = ctx->da2;
    PetscScalar     wx = 0.5 / ctx->hx, wy = 0.5 / ctx->hy;
    PetscInt        my = ctx->my;
    Vec             avec;
    PetscInt        i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ****a, ***b;

    DMGetLocalVector (da2, &avec);
    DMGlobalToLocalBegin (da2, Vvec, INSERT_VALUES, avec);
    DMGlobalToLocalEnd (da2, Vvec, INSERT_VALUES, avec);

    DMDAVecGetArrayDOFRead (da2, avec, &a);
    DMDAVecGetArray (da, bvec, &b);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            int         j0, j1;
            PetscScalar wyj;

            if (j == 0) {
                wyj = 2.0 * wy;
                j1  = j + 1;
                j0  = j; }
            else if (j == my - 1) {
                wyj = 2.0 * wy;
                j1  = j;
                j0  = j - 1; }
            else {
                wyj = wy;
                j1  = j + 1;
                j0  = j - 1; }

            for (i = xs; i < xs + xm; i++) {
                b[k][j][i] =
                    (a[k][j][i + 1][1] - a[k][j][i - 1][1]) * wx -
                    (a[k][j1][i][0] - a[k][j0][i][0]) * wyj; } } }

    DMDAVecRestoreArrayDOFRead (da, avec, &a);
    DMDAVecRestoreArray (da, bvec, &b);

    DMRestoreLocalVector (da, &avec);

    return (0); }


/* *
 * Horizontal advection operator, $ \mathbf{V} \cdot \nabla $
 *
 * - replaces the input field with the result
 * - one-sided derivatives at the north and south boundaries
 */

extern PetscErrorCode horizontal_advection (
    Vec bvec, Vec Vvec, Context ctx) {

    DM              da = ctx->da, da2 = ctx->da2;
    PetscScalar     wx = 0.5 / ctx->hx, wy = 0.5 / ctx->hy;
    PetscInt        my = ctx->my;
    Vec             avec;
    PetscInt        i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ****V, ***a, ***b;

    DMGetLocalVector (da, &avec);
    DMGlobalToLocalBegin (da, bvec, INSERT_VALUES, avec);
    DMGlobalToLocalEnd (da, bvec, INSERT_VALUES, avec);

    DMDAVecGetArrayRead (da, avec, &a);
    DMDAVecGetArrayDOFRead (da2, Vvec, &V);

    DMDAVecGetArray (da, bvec, &b);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            int         j0, j1;
            PetscScalar wyj;

            if (j == 0) {
                wyj = 2.0 * wy;
                j1  = j + 1;
                j0  = 0; }
            else if (j == my - 1) {
                wyj = 2.0 * wy;
                j1  = j;
                j0  = j - 1; }
            else {
                wyj = wy;
                j1  = j + 1;
                j0  = j - 1; }

            for (i = xs; i < xs + xm; i++) {
                b[k][j][i] =
                    V[k][j][i][0] * (a[k][j][i + 1] - a[k][j][i - 1]) *
                    wx +
                    V[k][j][i][1] * (a[k][j1][i] - a[k][j0][i]) * wyj; } } }

    DMDAVecRestoreArrayRead (da, avec, &a);
    DMDAVecRestoreArrayDOFRead (da, Vvec, &V);
    DMDAVecRestoreArray (da, bvec, &b);

    DMRestoreLocalVector (da, &avec);

    return (0); }


/* *
 * Operator coriolis parameter times pressure derivative,
 * $ f \frac{\partial}{\partial p} $
 *
 * - one-sided derivatives at top and bottom boundaries
 */

extern PetscErrorCode fpder (Vec bvec, Context ctx) {

    DM             da = ctx->da;
    PetscInt       mz = ctx->mz;
    PetscScalar    wz = 0.5 / ctx->hz;
    PetscScalar*   f  = ctx->Coriolis_parameter;
    Vec            avec;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ***a, ***b;

    DMGetLocalVector (da, &avec);
    DMGlobalToLocalBegin (da, bvec, INSERT_VALUES, avec);
    DMGlobalToLocalEnd (da, bvec, INSERT_VALUES, avec);


    DMDAVecGetArray (da, bvec, &b);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (k = zs; k < zs + zm; k++) {
        int         k0, k1;
        PetscScalar w;

        if (k == 0) {
            w  = 2.0 * wz;
            k1 = k + 1;
            k0 = 0; }
        else if (k == mz - 1) {
            w  = 2.0 * wz;
            k1 = k;
            k0 = k - 1; }
        else {
            w  = wz;
            k1 = k + 1;
            k0 = k - 1; }

        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                b[k][j][i] = f[j] * (a[k1][j][i] - a[k0][j][i]) * w; } } }

    DMDAVecRestoreArrayRead (da, avec, &a);
    DMDAVecRestoreArray (da, bvec, &b);

    DMRestoreLocalVector (da, &avec);

    return (0); }


extern PetscErrorCode horizontal_average (
    Context ctx, Vec v, PetscScalar v_ave[]) {

    DM             da;
    PetscInt       zs, ys, xs, zm, ym, xm, m, n;
    PetscScalar ***a, w;

    VecGetDM (v, &da);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    w = 1.0 / (double) ym / (double) xm;
    DMDAVecGetArrayRead (da, v, &a);

    for (int k   = 0; k < (int) ctx->mz; k++)
        v_ave[k] = 0.0;

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++)
                v_ave[k] += a[k][j][i] * w; } }

    DMDAVecRestoreArrayRead (da, v, &a);

#pragma push_macro("MPI_Allreduce")
#undef MPI_Allreduce
    MPI_Allreduce (
        MPI_IN_PLACE,
        v_ave,
        ctx->mz,
        MPI_DOUBLE,
        MPI_SUM,
        PETSC_COMM_WORLD);
#pragma pop_macro("MPI_Allreduce")

    DMDAGetInfo (da, 0, 0, 0, 0, &m, &n, 0, 0, 0, 0, 0, 0, 0);
    w = 1.0 / (double) m / (double) n;

    for (int k = 0; k < ctx->mz; k++)
        v_ave[k] *= w;

    return (0); }


/* *
 * R over pressure times laplace operator, $ \frac{R}{p} \nabla^2 $
 *
 * - one-sided derivatives at the north and south boundaries
 */

extern PetscErrorCode plaplace (Vec inout, Context ctx) {

    DM             da = ctx->da;
    PetscScalar*   p  = ctx->Pressure;
    PetscScalar    hx = ctx->hx;
    PetscScalar    hy = ctx->hy;
    PetscInt       my = ctx->my;
    PetscInt       mz = ctx->mz;
    const double   R  = Specific_gas_constant_of_dry_air;
    Vec            Vvec;
    PetscInt       zs, ys, xs, zm, ym, xm;
    PetscScalar*** v;
    PetscScalar*** result;

    DMGetLocalVector (da, &Vvec);
    DMGlobalToLocalBegin (da, inout, INSERT_VALUES, Vvec);
    DMGlobalToLocalEnd (da, inout, INSERT_VALUES, Vvec);

    DMDAVecGetArray (da, inout, &result);
    DMDAVecGetArray (da, Vvec, &v);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        double wx = R / (ctx->hx * ctx->hx * p[k]);
        double wy = R / (ctx->hy * ctx->hy * p[k]);

        for (int j = ys; j < ys + ym; j++) {
            int j0, j1;

            if (j == 0) {
                j1 = j;
                j0 = j; }
            else if (j == my - 1) {
                j1 = j;
                j0 = j; }
            else {
                j1 = j + 1;
                j0 = j - 1; }

            for (int i = xs; i < xs + xm; i++) {
                result[k][j][i] =
                    wx * (v[k][j][i + 1] - 2.0 * v[k][j][i] +
                          v[k][j][i - 1]) +
                    wy * (v[k][j1][i] - 2.0 * v[k][j][i] + v[k][j0][i]); } } }

    DMDAVecRestoreArray (da, inout, &result);
    DMDAVecRestoreArray (da, Vvec, &v);

    DMRestoreLocalVector (da, &Vvec);

    return (0); }
