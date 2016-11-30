#ifndef IO_H
#define IO_H

#include "petscdm.h"
#include "context.h"
#include "arrays.h"

extern const char *dimnames[4];

extern PetscErrorCode ncfile_open(const char *wrfin,int *ncid);

extern PetscErrorCode ncfile_get_dimsize(const int ncid,const char *dimname,
                                         size_t *dimsize);

/* Read 3d (lev,lat,lon) slice from 4d NetCDF array (Time,lev,lat,lon) */
extern PetscErrorCode read3D(const int ncid, const unsigned long time,
                             const char* varname, Vec v);

/* Read 2d (lat,lon) slice from 3d NetCDF array (Time,lat,lon) */
extern PetscErrorCode read2D(const int ncid, const unsigned long time,
                             const char* varname, Vec v);

/* Read 2d (lat,lon) slice from 4d NetCDF array (Time,lat,lon) */

extern PetscErrorCode readArray2D(const int ncid, const int time,
                                  const char *name,Array2D array);

extern PetscErrorCode readAttribute(
	const int    ncid,
	const char  *name,
	PetscScalar *attr);

extern PetscErrorCode readArray1D(
        const int            ncid,
        const unsigned long  time,
        const char          *varname,
        const int            n,
        PetscScalar         *a);

extern PetscErrorCode readArray(
        const int            ncid,
        const char          *varname,
        const size_t        *start,
        const size_t        *count,
        PetscScalar         *a);

extern PetscErrorCode write3D(const int ncid, const unsigned long time, \
                              const char *varname, Vec v);

extern PetscErrorCode write3Ddump(const char *varname, Vec v);

#endif
