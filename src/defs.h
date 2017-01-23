#ifndef UTILS_H
#define UTILS_H

#include <petscsys.h>
#include <petscvec.h>

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0])
#define ERROR(msg) SETERRQ (PETSC_COMM_WORLD, 1, msg)
#define WARNING(msg) PetscPrintf (PETSC_COMM_WORLD, "%s\n", msg)
#define DALOOP(LOOP_DA,LOOP_EXPRESSION)                       \
    do {                                                      \
        int LOOP_XS,LOOP_YS,LOOP_ZS,LOOP_XM,LOOP_YM,LOOP_ZM;  \
        DMDAGetCorners(LOOP_DA,&LOOP_XS,&LOOP_YS,             \
                       &LOOP_ZS,&LOOP_XM,&LOOP_YM,            \
                       &LOOP_ZM);                             \
        for (int k=LOOP_ZS; k<LOOP_ZS+LOOP_ZM; k++) {         \
            for (int j=LOOP_YS; j<LOOP_YS+LOOP_YM; j++) {     \
                for (int i=LOOP_XS; i<LOOP_XS+LOOP_XM;        \
                     i++) {LOOP_EXPRESSION;}}}} while(0)

#endif
