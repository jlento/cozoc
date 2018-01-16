#ifndef OPS_H
#define OPS_H

#include "context.h"
#include "petscvec.h"
#include <petscdmda.h>

int field_array1d_add (
    Vec x, PetscScalar* arr, DMDADirection direction);

int diff1d (
    const int n, PetscScalar* x, PetscScalar* f, PetscScalar* dfdx);

int horizontal_rotor (
    DM                da,
    DM                da2,
    const size_t      my,
    const PetscScalar hx,
    const PetscScalar hy,
    Vec               Vvec,
    Vec               bvec);

int horizontal_advection (Vec b, Vec V, Context* ctx);

int horizontal_average (Context* ctx, Vec v, PetscScalar v_ave[]);

int fpder (
    DM da, PetscInt mz, PetscScalar* f, PetscScalar* p, Vec bvec);

int plaplace (Vec inout, Context* ctx);
/*
int mul_fact (Context* ctx, Vec s);

int ellipticity_sigma_vorticity (
    Context*     ctx,
    size_t       mz,
    PetscScalar* p,
    PetscScalar* f,
    Vec          sigmavec,
    Vec          zetavec,
    Vec          V);

int xder (Vec bvec, Context* ctx);

int yder (Vec bvec, Context* ctx);
*/
#endif /* OPS_H */
