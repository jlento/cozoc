#ifndef OMEGA_H
#define OMEGA_H

#include "petscksp.h"

extern PetscErrorCode omega_compute_rhs_F_V(KSP ksp,Vec b,void *ctx_p);
extern PetscErrorCode omega_compute_operator(KSP ksp,Mat A,Mat B,void *ctx_p);


#endif /* OMEGA_H */
