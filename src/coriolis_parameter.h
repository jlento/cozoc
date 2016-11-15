#ifndef CORIOLIS_PARAMETER_H
#define CORIOLIS_PARAMETER_H

#include "petscvec.h"
#include "context.h"

extern PetscErrorCode coriolis_parameter_read(Context ctx);
extern PetscErrorCode coriolis_parameter_add(Vec x,Context ctx);


#endif /* CORIOLIS_PARAMETER_H */
