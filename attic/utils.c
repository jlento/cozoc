#include "utils.h"
#include <petscerror.h>

int norm2_and_absmax (const char* name, Vec x) {
    PetscScalar norm;

    VecNorm (x, NORM_2, &norm);
    PetscPrintf (
        PETSC_COMM_WORLD,
        "%s's norm and |max|: %12g  ",
        name,
        (double) norm);
    VecNorm (x, NORM_INFINITY, &norm);
    PetscPrintf (PETSC_COMM_WORLD, "%12g\n", (double) norm);

    return (0); }
