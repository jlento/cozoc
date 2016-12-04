#ifndef UTILS_H
#define UTILS_H

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0])

#include <petscsys.h>
#include <petscvec.h>
extern PetscErrorCode norm2_and_absmax(const char *name,Vec x);

#endif

