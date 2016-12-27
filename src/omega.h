#ifndef OMEGA_H
#define OMEGA_H

#include "petscksp.h"

enum { V, T, F, Q, N_OMEGA_COMPONENTS};

extern char* omega_component_id_string[N_OMEGA_COMPONENTS];

extern PetscErrorCode omega_compute_operator(KSP ksp,Mat A,Mat B,void *ctx_p);

extern PetscErrorCode omega_compute_rhs_F_V(KSP ksp,Vec b,void *ctx_p);
extern PetscErrorCode omega_compute_rhs_F_T(KSP ksp,Vec b,void *ctx_p);
extern PetscErrorCode omega_compute_rhs_F_F(KSP ksp,Vec b,void *ctx_p);
extern PetscErrorCode omega_compute_rhs_F_Q(KSP ksp,Vec b,void *ctx_p);

extern PetscErrorCode (*omega_compute_rhs[N_OMEGA_COMPONENTS])(KSP ksp,Vec b,void *ctx_p);



#endif /* OMEGA_H */
