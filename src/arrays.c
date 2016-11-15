#include "arrays.h"
#include "petscsys.h"

#undef __FUNCT__
#define __FUNCT__ "allocateArray2D"
extern PetscErrorCode allocateArray2D(const int xs,const int ys,
                                      const int xm,const int ym,
                                      Array2D *array)
{
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        array->xs = xs; array->xm = xm;
        array->ys = ys; array->ym = ym;
        //array->data  = malloc(array->xm*array->ym*sizeof(PetscScalar));
        ierr = PetscMalloc1(array->xm*array->ym,&array->data);CHKERRQ(ierr);
        PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "freeArray2D"
                extern PetscErrorCode freeArray2D(Array2D *array)
                {
                        PetscFunctionBeginUser;
                        free(array->data);
                        PetscFunctionReturn(0);
                }

#undef __FUNCT__
#define __FUNCT__ "arrayGetArray2D"
        extern PetscErrorCode arrayGetArray2D(Array2D array, PetscScalar **a[])
        {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;
                ierr = PetscMalloc1(array.ym,a);CHKERRQ(ierr);
                for (int i=0; i<array.ym; i++)
                        (*a)[i] = array.data + i*array.xm - array.xs;
                *a -= array.ys;
                PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "arrayRestoreArray2D"
        extern PetscErrorCode arrayRestoreArray2D(Array2D array, PetscScalar **a[])
        {
                void *tmp;
                PetscErrorCode ierr;
                PetscFunctionBeginUser;
                tmp = (void*)(*a + array.ys);
                ierr = PetscFree(tmp);CHKERRQ(ierr);
                PetscFunctionReturn(0);
        }
