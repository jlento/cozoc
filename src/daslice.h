#ifndef DASLICE_H
#define DASLICE_H

#include "petscdmda.h"

extern PetscErrorCode create_subdm_plane(DMDADirection dir,DM da,DM *dap);


#endif /* DASLICE_H */
