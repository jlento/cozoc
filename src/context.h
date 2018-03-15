#pragma once

#include "io.h"
#include "omega.h"
#include "options.h"
#include <petscksp.h>

typedef struct Context Context;
struct Context {
    size_t       current, first, last;
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    PetscInt     mx, my, mz, mt;    // Global grid sizes
    PetscScalar  hx, hy, hz;        // Grid spacings
    int          cu_physics;
    double *     Time_coordinate;    // In seconds
    PetscScalar *Pressure;
    PetscScalar *Coriolis_parameter;
    Vec          Diabatic_heating;
    Vec          Friction;
    Vec          Geopotential_height;
    Vec          Horizontal_wind;
    Vec          omega[NUM_GENERALIZED_OMEGA_COMPONENTS];
    Vec          One_over_dry_air_mass_column;
    Vec          Temperature;
    Vec          Temperature_tendency;
    Vec          Sigma_parameter;
    Vec          Surface_pressure;
    Vec          Vorticity;
    Vec          Vorticity_tendency;
};

Context new_context (Options, NCFile);

Vec new_vec (Context *);

void free_context (Context *ctx);

void update_context (size_t, NCFile, Context *);

int diabatic_heating (Context *, const int, const int);
