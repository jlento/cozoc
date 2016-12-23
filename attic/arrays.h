#ifndef ARRAYS_H
#define ARRAYS_H

#include <petscsys.h>

typedef struct {
        int xs,ys,xm,ym;
        PetscScalar *data;
} Array2D;

extern PetscErrorCode allocateArray2D(const int xs,const int ys,
                                      const int xm,const int ym,
                                      Array2D *array);
extern PetscErrorCode freeArray2D(Array2D *array);
extern PetscErrorCode arrayGetArray2D(Array2D array, PetscScalar ***a);
extern PetscErrorCode arrayRestoreArray2D(Array2D array, PetscScalar ***a);

#endif /* ARRAYS_H */
