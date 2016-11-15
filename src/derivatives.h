#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "petscvec.h"
#include "context.h"
extern PetscErrorCode pder(const PetscScalar hz,const Vec fvec,Vec dvec);
extern PetscErrorCode horizontal_rotor(Vec Vvec,Vec bvec,Context ctx);
extern PetscErrorCode horizontal_advection(Vec V,Vec b,Context ctx);
extern PetscErrorCode fpder(Vec b,Context ctx);

#endif /* DERIVATIVES_H */
