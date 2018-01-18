#ifndef CONTEXT_H
#define CONTEXT_H

#include "omega.h"
#include "options.h"
#include <petscksp.h>

typedef struct {
    int ncid;
    int skip;
    int steps;
    int flags;
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    size_t       mx, my, mz, mt;    // Global grid sizes
    PetscScalar  hx, hy, hz;        // Grid spacings
    int          cu_physics;
    double*      Time_coordinate; // In seconds
    PetscScalar* Pressure;
    PetscScalar* Coriolis_parameter;
    Vec          Surface_pressure;
    Vec          Temperature;
    Vec          Sigma_parameter;
    Vec          Vorticity;
    Vec          Geopotential_height;
    Vec          Diabatic_heating;
    Vec          Horizontal_wind;
    Vec          Friction;
    Vec          Temperature_tendency;
    Vec          Vorticity_tendency;
    Vec          omega[N_OMEGA_COMPONENTS]; } Context;


int context_create (
    const int ncid, PetscInt, PetscInt, Context*);

int context_destroy (Context* ctx);

int context_update (const int ncid, const int step, Context* ctx);

#endif /* CONTEXT_H */
