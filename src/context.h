#ifndef CONTEXT_H
#define CONTEXT_H

#include "io.h"
#include "omega.h"
#include "options.h"
#include <petscksp.h>

typedef struct nContext nContext;
struct nContext {
    size_t       step;
    size_t       dim[NUM_DIM];
    PetscScalar  delta[NUM_DIM];
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    size_t       mx, my, mz, mt;    // Global grid sizes
    PetscScalar  hx, hy, hz;        // Grid spacings
    int          cu_physics;
    double *     Time_coordinate;    // In seconds
    PetscScalar *Pressure;
    PetscScalar *Coriolis_parameter;
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
    Vec          omega[NUM_GENERALIZED_OMEGA_COMPONENTS];
};

nContext new_context (Options, NCFile);

typedef struct {
    int          ncid;
    size_t       first;
    size_t       last;
    int          flags;
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    size_t       mx, my, mz, mt;    // Global grid sizes
    PetscScalar  hx, hy, hz;        // Grid spacings
    int          cu_physics;
    double *     Time_coordinate;    // In seconds
    PetscScalar *Pressure;
    PetscScalar *Coriolis_parameter;
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
    Vec          omega[NUM_GENERALIZED_OMEGA_COMPONENTS];
} Context;

int context_create (NCFile, Context *);

int context_destroy (Context *ctx);

int context_update (
    NCFile, size_t const step, size_t const first, size_t const last,
    Context *ctx);

#endif /* CONTEXT_H */
