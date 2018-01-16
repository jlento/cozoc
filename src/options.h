#pragma once

#include <petscsys.h>


typedef struct Options Options;
struct Options {
    char fname[PETSC_MAX_PATH_LEN];
    PetscInt skip;
    PetscInt steps;
    PetscBool compute_omega_quasi_geostrophic;
    PetscBool compute_omega_generalized;
};

extern PetscErrorCode read_options (Options*);
