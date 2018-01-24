#ifndef OMEGA_H
#define OMEGA_H

#include "petscksp.h"
#include <netcdf.h>

typedef enum GENERALIZED_OMEGA_COMPONENTS {
    GENERALIZED_OMEGA_COMPONENT_V,
    GENERALIZED_OMEGA_COMPONENT_T,
    GENERALIZED_OMEGA_COMPONENT_F,
    GENERALIZED_OMEGA_COMPONENT_Q,
    GENERALIZED_OMEGA_COMPONENT_A
} GENERALIZED_OMEGA_COMPONENTS;

#define NUM_GENERALIZED_OMEGA_COMPONENTS 5

extern char omega_component_id_string[NUM_GENERALIZED_OMEGA_COMPONENTS]
                                     [NC_MAX_NAME];

extern PetscErrorCode omega_compute_operator (KSP ksp, Mat A, Mat B, void *ctx);

extern PetscErrorCode omega_compute_rhs_F_V (KSP ksp, Vec b, void *ctx);
extern PetscErrorCode omega_compute_rhs_F_T (KSP ksp, Vec b, void *ctx);
extern PetscErrorCode omega_compute_rhs_F_F (KSP ksp, Vec b, void *ctx);
extern PetscErrorCode omega_compute_rhs_F_Q (KSP ksp, Vec b, void *ctx);
extern PetscErrorCode omega_compute_rhs_F_A (KSP ksp, Vec b, void *ctx);

extern PetscErrorCode (*omega_compute_rhs[NUM_GENERALIZED_OMEGA_COMPONENTS]) (
    KSP ksp, Vec b, void *ctx);

#endif /* OMEGA_H */
