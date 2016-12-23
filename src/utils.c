#include "utils.h"
#include <petscerror.h>

#undef __FUNCT__
#define __FUNCT__ "norm2_and_absmax"
extern PetscErrorCode norm2_and_absmax(const char *name, Vec x) {
    PetscScalar    norm;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecNorm(x, NORM_2, &norm);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s's norm and |max|: %12g  ", name,
                       (double)norm);
    CHKERRQ(ierr);
    ierr = VecNorm(x, NORM_INFINITY, &norm);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%12g\n", (double)norm);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
