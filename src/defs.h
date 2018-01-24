#pragma once

#include <petscsys.h>
#include <petscdmda.h>

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0])
#define ERROR(msg) SETERRQ (PETSC_COMM_WORLD, 1, msg)
#define WARNING(msg) PetscPrintf (PETSC_COMM_WORLD, "%s\n", msg)

/*
typedef enum FIELD FIELD;
enum FIELD { TEMPERATURE };

typedef Vec (*UpdateFunction)(FIELD, size_t, int);

typedef struct Update Update;
struct Update {
    UpdateFunction temperature;
};

typedef Update (*RegisterUpdateFunctionsFunction)(void);
*/
