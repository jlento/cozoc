#include "vecops.h"
#include "context.h"
#include "io.h"
#include "petscdmda.h"
#include <strings.h>



#undef __FUNCT__
#define __FUNCT__ "xder"
extern PetscErrorCode xder(Vec f, const PetscScalar dx, Vec dfdx) {
    DM             da;
    PetscInt       i, j, k, xs, ys, zs, xm, ym, zm;
    Vec            floc;
    PetscScalar ***fa, ***dfdxa;
    PetscScalar    w = 0.5 / dx;
    PetscErrorCode ierr;

    ierr = VecGetDM(f, &da);
    CHKERRQ(ierr);

    ierr = DMGetLocalVector(da, &floc);
    CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da, f, INSERT_VALUES, floc);
    CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, f, INSERT_VALUES, floc);
    CHKERRQ(ierr);

    ierr = DMDAVecGetArray(da, floc, &fa);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, dfdx, &dfdxa);
    CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    CHKERRQ(ierr);
    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                dfdxa[k][j][i] = (fa[k][j][i + 1] - fa[k][j][i - 1]) * w; } } }
    ierr = DMDAVecRestoreArray(da, floc, &fa);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, dfdx, &dfdxa);
    CHKERRQ(ierr);

    ierr = DMRestoreLocalVector(da, &floc);
    CHKERRQ(ierr);

    PetscFunctionReturn(0); }

#undef __FUNCT__
#define __FUNCT__ "yder"
extern PetscErrorCode yder(Vec f, const PetscScalar dy, Vec dfdy) {
    DM             da;
    Vec            floc;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm, my;
    PetscScalar    w;
    PetscScalar ***fa, ***dfdya;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecGetDM(f, &da);
    CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    CHKERRQ(ierr);
    ierr = DMDAGetInfo(da, 0, 0, &my, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    CHKERRQ(ierr);

    ierr = DMCreateLocalVector(da, &floc);
    CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da, f, INSERT_VALUES, floc);
    CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, f, INSERT_VALUES, floc);
    CHKERRQ(ierr);

    ierr = DMDAVecGetArray(da, floc, &fa);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, dfdy, &dfdya);
    CHKERRQ(ierr);
    w = 1.0 / dy;
    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                if (j == 0) {
                    dfdya[k][j][i] = (fa[k][j + 1][i] - fa[k][j][i]) * w; }
                else if (j == my - 1) {
                    dfdya[k][j][i] = (fa[k][j][i] - fa[k][j - 1][i]) * w; }
                else {
                    dfdya[k][j][i] =
                        (fa[k][j + 1][i] - fa[k][j - 1][i]) * w * 0.5; } } } }
    ierr = DMDAVecRestoreArray(da, dfdy, &dfdya);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, floc, &fa);
    CHKERRQ(ierr);

    ierr = VecDestroy(&floc);
    CHKERRQ(ierr);

    PetscFunctionReturn(0); }

#undef __FUNCT__
#define __FUNCT__ "curl"
extern PetscErrorCode curl(Vec u, Vec v, PetscScalar dx, PetscScalar dy,
                           Vec f) {
    Vec            dudy, dvdx;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    dvdx = f;
    ierr = VecDuplicate(dvdx, &dudy);
    CHKERRQ(ierr);

    ierr = yder(u, dy, dudy);
    CHKERRQ(ierr);
    ierr = xder(v, dx, dvdx);
    CHKERRQ(ierr);
    ierr = VecAXPY(dvdx, -1.0, dudy);
    CHKERRQ(ierr);

    ierr = VecDestroy(&dudy);
    CHKERRQ(ierr);

    PetscFunctionReturn(0); }

#undef __FUNCT__
#define __FUNCT__ "Vec2Vec123AXPY"
extern PetscErrorCode Vec2Vec123AXPY(Vec y, PetscScalar a, Vec x) {
    DM             da;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ***ya, *xa;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecGetDM(y, &da);
    CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    CHKERRQ(ierr);

    ierr = DMDAVecGetArray(da, y, &ya);
    CHKERRQ(ierr);
    ierr = VecGetArray(x, &xa);
    CHKERRQ(ierr);

    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                ya[k][j][i] += a * xa[j]; } } }

    ierr = VecRestoreArray(x, &xa);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, y, &ya);
    CHKERRQ(ierr);

    PetscFunctionReturn(0); }

#undef __FUNCT__
#define __FUNCT__ "advection"
extern PetscErrorCode advection(Vec u, Vec v, Vec f, PetscScalar dx,
                                PetscScalar dy, Vec adv) {
    Vec            fx, fy;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    fy   = adv;
    ierr = VecDuplicate(fy, &fx);
    CHKERRQ(ierr);

    /* fx = df/dx; fy = df/dy */
    ierr = xder(f, dx, fx);
    CHKERRQ(ierr);
    ierr = yder(f, dy, fy);
    CHKERRQ(ierr);

    /* fx = u * df/dx; fy = v * df/dy */
    ierr = VecPointwiseMult(fx, u, fx);
    CHKERRQ(ierr);
    ierr = VecPointwiseMult(fy, v, fy);
    CHKERRQ(ierr);

    /* adv = u * df/dx + v * df/dy */
    ierr = VecAXPY(fy, 1.0, fx);
    CHKERRQ(ierr);

    ierr = VecDestroy(&fx);
    CHKERRQ(ierr);

    PetscFunctionReturn(0); }

#undef __FUNCT__
#define __FUNCT__ "ydivide"
extern PetscErrorCode ydivide(Vec w, Vec x, Vec y) {
    DM             da;
    PetscInt       i, j, k, zs, ys, xs, zm, ym, xm;
    PetscScalar ***wa, ***xa, *ya;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecGetDM(w, &da);
    CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, w, &wa);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, x, &xa);
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &ya);
    CHKERRQ(ierr);
    for (k = zs; k < zs + zm; k++) {
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                wa[k][j][i] = xa[k][j][i] / ya[j]; } } }
    ierr = DMDAVecRestoreArray(da, w, &wa);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, x, &xa);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &ya);
    CHKERRQ(ierr);
    PetscFunctionReturn(0); }
