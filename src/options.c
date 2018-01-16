#include "options.h"
#include <petscoptions.h>

PetscErrorCode read_options (Options *options) {
    PetscOptionsBegin (PETSC_COMM_WORLD, "", "Options for COZOC", "none");

    CHKERRQ (PetscOptionsString ("-f",
        "Input file, NetCDF4/HDF5 format, from WRF simulation", 0,
        options->fname, options->fname, PETSC_MAX_PATH_LEN, 0));

    CHKERRQ (PetscOptionsInt (
        "-s", "Skip <n> timesteps", 0, options->skip, &options->skip, 0));
    CHKERRQ (PetscOptionsInt ("-n", "Calculate <n> timesteps", 0,
        options->steps, &options->steps, 0));

    CHKERRQ (PetscOptionsBool ("-Q", "Calculate quasi-geostrophic omega eq.", 0,
        options->compute_omega_quasi_geostrophic,
        &options->compute_omega_quasi_geostrophic, 0));
    CHKERRQ (PetscOptionsBool ("-G", "Calculate generalized omega eq.", 0,
        options->compute_omega_generalized, &options->compute_omega_generalized,
        0));
    PetscOptionsEnd ();
    return 0;
};
