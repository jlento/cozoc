#ifndef UTILS_H
#define UTILS_H

#include <petscsys.h>

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0])
#define ERROR(msg) SETERRQ (PETSC_COMM_WORLD, 1, msg)
#define WARNING(msg) PetscPrintf (PETSC_COMM_WORLD, "%s\n", msg)

#endif
