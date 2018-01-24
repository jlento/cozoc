#ifndef OMEGAQG_H
#define OMEGAQG_H

#include "petscksp.h"
#include "context.h"

#define OMEGA_QG_ID_STRING "ome_qg"

#define NUM_QUASI_GEOSTROPHIC_OMEGA_EQS 1

extern PetscErrorCode omega_qg_compute_operator(KSP ksp, Mat dummy,
                                                Mat L, void *ctx_p);

extern PetscErrorCode omega_qg_compute_rhs(KSP ksp,Vec b,void *ctx_p);

#endif /* OMEGAQG_H */
