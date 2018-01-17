#include "options.h"
#include "abortonerror.h"
#include <petscoptions.h>

Options new_options () {

    Options options = {.fname = "wrf.nc4",
                       .skip = 0,
                       .steps = PETSC_MAX_INT,
                       .compute_omega_quasi_geostrophic = PETSC_TRUE,
                       .compute_omega_generalized = PETSC_TRUE};

    PetscOptionsBegin (PETSC_COMM_WORLD, "", "Options for COZOC", "none");

    CHKERRQ (
        PetscOptionsString (
            "-f", "Input file, NetCDF4/HDF5 format, from WRF simulation", 0,
            options.fname, options.fname, PETSC_MAX_PATH_LEN, 0));

    CHKERRQ (
        PetscOptionsInt (
            "-s", "Skip <n> timesteps", 0, options.skip, &options.skip, 0));
    CHKERRQ (
        PetscOptionsInt (
            "-n", "Calculate <n> timesteps", 0, options.steps, &options.steps,
            0));

    CHKERRQ (
        PetscOptionsBool (
            "-Q", "Disable quasi-geostrophic omega eq. calculation", 0,
            options.compute_omega_quasi_geostrophic,
            &options.compute_omega_quasi_geostrophic, 0));
    CHKERRQ (
        PetscOptionsBool (
            "-G", "Disable generalized omega eq. calculation", 0,
            options.compute_omega_generalized,
            &options.compute_omega_generalized, 0));

    PetscOptionsEnd ();

    return options;
};
