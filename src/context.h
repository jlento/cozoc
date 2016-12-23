#ifndef CONTEXT_H
#define CONTEXT_H

#include "omega.h"
#include <petscksp.h>

enum {
    OMEGA_QUASI_GEOSTROPHIC = (1 << 0),
    OMEGA_GENERALIZED       = (1 << 1),
    HEIGHT_TENDENCY         = (1 << 2) } Targets;

typedef struct {
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    size_t       mx, my, mz, mt; /* Global grid sizes */
    PetscScalar  hx, hy, hz;     /* Grid spacings */
    PetscInt     time;           /* Current timestep */
    PetscScalar* Pressure;
    PetscScalar* Coriolis_parameter;
    Vec          Surface_pressure;
    Vec          Sigma_parameter;
    Vec          Vorticity;
    Vec          Horizontal_wind;
    Vec          Geopotential_height;
    Vec          Temperature;
    Vec          omega[N_OMEGA_COMPONENTS]; } tContext, *Context;


extern PetscErrorCode context_create (
    const int ncid, int skip, int* steps, int* flags, Context* ctx);

extern PetscErrorCode context_destroy (Context* ctx);

extern PetscErrorCode context_update (
    const int ncid, const int time, Context ctx);

#endif /* CONTEXT_H */
