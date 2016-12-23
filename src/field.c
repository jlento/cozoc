#include "field.h"
#include <petscdmda.h>
#include <stdarg.h>


extern PetscErrorCode field_array1d_add (
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


extern PetscErrorCode field_array2d_add (
    Vec x, PetscScalar** arr, DMDADirection direction) {

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
                    xa[k][j][i] += arr[j][k]; } } }

    case DMDA_Y:
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    xa[k][j][i] += arr[i][k]; } } }

    case DMDA_Z:
        for (k = zs; k < zs + zm; k++) {
            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    xa[k][j][i] += arr[i][j]; } } } }

    DMDAVecRestoreArray (da, x, &xa);
    return (0); }
