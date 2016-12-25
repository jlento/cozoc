#ifndef UTILS_H
#define UTILS_H

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0])
#define ERROR(msg) SETERRQ (PETSC_COMM_WORLD, 1, msg)

#include <petscsys.h>
#include <petscvec.h>

int norm2_and_absmax (const char* name, Vec x);

#endif
