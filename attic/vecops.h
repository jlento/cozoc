#ifndef VECOPS_H
#define VECOPS_H

#include "petscvec.h"
#include "context.h"

extern PetscErrorCode horizontal_average(Context,Vec,PetscScalar*);
extern PetscErrorCode xder(Vec f,const PetscScalar dx,Vec dfdx);
extern PetscErrorCode yder(Vec f,const PetscScalar dy,Vec dfdy);
extern PetscErrorCode curl(Vec u,Vec v,PetscScalar dx,PetscScalar dy,Vec f);
extern PetscErrorCode Vec2Vec123AXPY(Vec y,PetscScalar a,Vec x);
extern PetscErrorCode advection(Vec u,Vec v,Vec f,
                                PetscScalar dx,PetscScalar dy,Vec adv);
extern PetscErrorCode ydivide(Vec w,Vec x,Vec y);

#endif /* VECOPS_H */
