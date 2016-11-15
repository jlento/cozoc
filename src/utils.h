#ifndef UTILS_H
#define UTILS_H

#include <petscsys.h>
#include <petscvec.h>

/* NetCDF error handler */
#define ERRCODE 2
#define ERR(e) if (e) {printf("%s[%d]: %s: Error: %s\n", __FILE__, __LINE__, __func__, nc_strerror(e)); exit(ERRCODE);}

extern PetscErrorCode norm2_and_absmax(const char *name,Vec x);

#endif
