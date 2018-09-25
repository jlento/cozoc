#pragma once

#include "io.h"
#include "omega.h"
#include "options.h"
#include <petscksp.h>

typedef struct Context Context;
struct Context {
    int          ncid; // Input file
    size_t       current, first, last;
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    PetscInt     mx, my, mz, mt;    // Global grid sizes
    PetscScalar  hx, hy, hz;        // Grid spacings
    double *     Time_coordinate;    // In seconds
    PetscScalar *Pressure;
    PetscScalar *Coriolis_parameter;
    Vec          Diabatic_heating;
    Vec          Diabatic_heating_attennuated;
    Vec          Diabatic_heating_forcing;
    Vec          Friction;
    Vec          Geopotential_height;
    Vec          Horizontal_wind;
    Vec          omega[NUM_GENERALIZED_OMEGA_COMPONENTS];
    Vec          One_over_dry_air_mass_column;
    Vec          Temperature;
    Vec          Temperature_tendency;
    Vec          Sigma_parameter;
    Vec          Surface_pressure;
    Vec          Surface_attennuation;
    Vec          Vorticity;
    Vec          Vorticity_tendency;
};

Context new_context (Options, Files);

Vec new_vec (Context *);

void free_context (Context *ctx);

void update_context (size_t, Files, Context *);

int diabatic_heating (Context *, const int ncid, const int time);
int friction (Context *, const int ncid, const int time);
int horizontal_wind_and_vorticity_and_vorticity_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, DM da, DM da2,
    size_t my, PetscScalar hx, PetscScalar hy, Vec V, Vec zeta,
    Vec zetatend, Context *ctx);
int one_over_dry_air_mass_column (const int ncid, const int time, Context *);
int temperature (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec T,
    Vec Ttend, Context *);
int sigma_parameter (
    DM da, PetscInt mz, PetscScalar *p, Vec Tvec, Vec sigmavec);
