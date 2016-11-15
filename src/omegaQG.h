#ifndef OMEGAQG_H
#define OMEGAQG_H

#include "petscksp.h"
#include "context.h"

extern PetscErrorCode compute_operator_omega_QG(KSP ksp, Mat dummy,
                                                Mat L, void *ctx_p);

extern PetscErrorCode compute_rhs_omega_QG(KSP ksp,Vec b,void *ctx_p);

#endif /* OMEGAQG_H */
