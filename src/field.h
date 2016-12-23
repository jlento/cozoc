#ifndef FIELD_H
#define FIELD_H

#include <petscdmda.h>

#define FIELD_ARRAY_ADD(x, arr, direction)                              \
    _Generic((arr),                                                     \
             PetscScalar*: field_array1d_add,                           \
             PetscScalar** : field_array2d_add) (x, arr, direction)

extern PetscErrorCode field_array1d_add (
    Vec x, PetscScalar* arr, DMDADirection direction);

extern PetscErrorCode field_array2d_add (
    Vec x, PetscScalar** arr, DMDADirection direction);

#endif /* FIELD_H */
