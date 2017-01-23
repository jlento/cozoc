#include "constants.h"
#include "context.h"
#include "ops.h"
#include "petscdmda.h"


int field_array1d_add (
    Vec x, PetscScalar* arr, DMDADirection direction) {

    DM             da;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar*** xa;

    VecGetDM (x, &da);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);
    DMDAVecGetArray (da, x, &xa);

    switch (direction) {
    case DMDA_X:
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    xa[k][j][i] += arr[i]; } } }

    case DMDA_Y:
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    xa[k][j][i] += arr[j]; } } }

    case DMDA_Z:
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    xa[k][j][i] += arr[k]; } } } }

    DMDAVecRestoreArray (da, x, &xa);
    return (0); }


int diff1d (
    const int n, PetscScalar* x, PetscScalar* f, PetscScalar* dfdx) {

    dfdx[0] = (f[1] - f[0]) / (x[1] - x[0]);

    for (int i  = 1; i < n - 1; i++)
        dfdx[i] = (f[i + 1] - f[i - 1]) / (x[i + 1] - x[i - 1]);

    dfdx[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);

    return (0); }

/* *
 * Horizontal rotor operator, $ \nabla \times $
 *
 * - one-sided derivatives at the north and south boundaries
 */

int horizontal_rotor (
    DM                da,
    DM                da2,
    const size_t      my,
    const PetscScalar hx,
    const PetscScalar hy,
    Vec               Vvec,
    Vec               bvec) {

    PetscScalar     wx = 0.5 / hx, wy = 0.5 / hy;
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
            else if (j == (int)my - 1) {
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

int horizontal_advection (Vec bvec, Vec Vvec, Context* ctx) {

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
 * Operator coriolis parameter (or some other f[y]) times pressure
 * derivative,
 * $ f \frac{\partial}{\partial p} $
 *
 * - one-sided derivatives at top and bottom boundaries
 */

int fpder (
    DM da, PetscInt mz, PetscScalar* f, PetscScalar* p, Vec bvec) {

    PetscScalar    wz = 0.5 / (p[1] - p[0]);
    Vec            avec;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ***a, ***b;

    DMGetLocalVector (da, &avec);
    DMGlobalToLocalBegin (da, bvec, INSERT_VALUES, avec);
    DMGlobalToLocalEnd (da, bvec, INSERT_VALUES, avec);

    DMDAVecGetArrayRead (da, avec, &a);
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

        if (f == NULL)
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    b[k][j][i] = (a[k1][j][i] - a[k0][j][i]) * w; } }
        else
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    b[k][j][i] = f[j] * (a[k1][j][i] - a[k0][j][i]) * w; } } }

    DMDAVecRestoreArrayRead (da, avec, &a);
    DMDAVecRestoreArray (da, bvec, &b);

    DMRestoreLocalVector (da, &avec);

    return (0); }


int horizontal_average (Context* ctx, Vec v, PetscScalar v_ave[]) {

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

    for (int k = 0; k < (int)ctx->mz; k++)
        v_ave[k] *= w;

    return (0); }


/* *
 * R over pressure times laplace operator, $ \frac{R}{p} \nabla^2 $
 *
 * - one-sided derivatives at the north and south boundaries
 */

int plaplace (Vec inout, Context* ctx) {

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
        double wx = R / (hx * hx * p[k]);
        double wy = R / (hy * hy * p[k]);

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
