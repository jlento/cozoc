#pragma once

#include <petscsys.h>

typedef struct Options Options;
struct Options {
    char      infname[PETSC_MAX_PATH_LEN];
    char      outfname[PETSC_MAX_PATH_LEN];
    size_t    first;
    size_t    last;
    PetscBool compute_omega_quasi_geostrophic;
    PetscBool compute_omega_generalized;
};

Options new_options (void);
