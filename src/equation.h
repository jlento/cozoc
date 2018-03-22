#pragma once

#include "omega.h"
#include "omegaQG.h"
#include <netcdf.h>
#include <petscksp.h>

typedef PetscErrorCode ComputeLHSOperatorFunction (KSP, Mat, Mat, void *ctx);
typedef PetscErrorCode ComputeRHSVectorFunction (KSP, Vec, void *ctx);

typedef struct Equations Equations;
struct Equations {
    size_t num_eq;
    ComputeLHSOperatorFunction
        *L[NUM_QUASI_GEOSTROPHIC_OMEGA_EQS + NUM_GENERALIZED_OMEGA_COMPONENTS];
    ComputeRHSVectorFunction
        *a[NUM_QUASI_GEOSTROPHIC_OMEGA_EQS + NUM_GENERALIZED_OMEGA_COMPONENTS];
    char id_string[NUM_QUASI_GEOSTROPHIC_OMEGA_EQS +
                   NUM_GENERALIZED_OMEGA_COMPONENTS][NC_MAX_NAME];
};

Equations new_equations (Options);
Vec solution (ComputeLHSOperatorFunction, ComputeRHSVectorFunction, Context);
