#include "options.h"
#include "abortonerror.h"
#include <petscoptions.h>

#define MAXLEN 256

Options new_options () {

    Options options = {.fname = "wrf.nc4",
                       .first = 0,
                       .last = PETSC_MAX_INT,
                       .compute_omega_quasi_geostrophic = PETSC_TRUE,
                       .compute_omega_generalized = PETSC_TRUE};
    char s[MAXLEN] = "";

    PetscOptionsBegin (PETSC_COMM_WORLD, "", "Options for COZOC", "none");

    CHKERRQ (
        PetscOptionsString (
            "-f", "Input file, NetCDF4/HDF5 format, from WRF simulation", 0,
            options.fname, options.fname, PETSC_MAX_PATH_LEN, 0));

    CHKERRQ (
        PetscOptionsString (
            "-r",
            "Range of steps to compute, counting form zero, -r <start>,<stop>",
            0, s, s, MAXLEN, 0));
    sscanf (s, "%u,%u", &options.first, &options.last);

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
