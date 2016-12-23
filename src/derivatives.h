#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "context.h"
#include "petscvec.h"

extern PetscErrorCode diff1d (
    const int n, PetscScalar* x, PetscScalar* f, PetscScalar* dfdx);

extern PetscErrorCode pder (
    const PetscScalar hz, const Vec fvec, Vec dvec);

extern PetscErrorCode horizontal_rotor (
    Vec Vvec, Vec bvec, Context ctx);

extern PetscErrorCode horizontal_advection (Vec b, Vec V, Context ctx);

extern PetscErrorCode horizontal_average (
    Context ctx, Vec v, PetscScalar v_ave[]);

extern PetscErrorCode fpder (Vec b, Context ctx);

extern PetscErrorCode plaplace (Vec inout, Context ctx);

#endif /* DERIVATIVES_H */
