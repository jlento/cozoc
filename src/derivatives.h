#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "context.h"
#include "petscvec.h"

int diff1d (
    const int n, PetscScalar* x, PetscScalar* f, PetscScalar* dfdx);

int pder (const PetscScalar hz, const Vec fvec, Vec dvec);

int horizontal_rotor (Vec Vvec, Vec bvec, Context ctx);

int horizontal_advection (Vec b, Vec V, Context ctx);

int horizontal_average (Context ctx, Vec v, PetscScalar v_ave[]);

int fpder (Vec b, Context ctx);

int plaplace (Vec inout, Context ctx);

#endif /* DERIVATIVES_H */
