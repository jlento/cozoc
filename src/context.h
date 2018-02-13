#pragma once

#include "omega.h"
#include "options.h"
#include "io.h"
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
    Vec          diab;
    Vec          rthcuten;
    Vec          rthraten;
    Vec          rthblten;
};

Context new_context (Options, NCFile);

Vec new_vec(Context*);

void free_context (Context *ctx);

void update_context (size_t, NCFile, Context *);
