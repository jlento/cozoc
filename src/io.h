#ifndef IO_H
#define IO_H

#include "context.h"
#include "netcdf.h"
#include "petscdm.h"


/* NetCDF error handler */
#define ERRCODE 2
#define ERR(e)                                                         \
    if (e) {                                                           \
        printf (                                                       \
            "%s[%d]: %s: Error: %s\n",                                 \
            __FILE__,                                                  \
            __LINE__,                                                  \
            __func__,                                                  \
            nc_strerror (e));                                          \
        exit (ERRCODE);                                                \
    }


extern PetscErrorCode file_open (const char* wrfin, int* ncid);
extern PetscErrorCode file_close (const int ncid);
extern PetscErrorCode file_redef (const int ncid);
extern PetscErrorCode file_enddef (const int ncid);
extern PetscErrorCode file_def_var (const int ncid, const char* name);
extern PetscErrorCode file_get_dimsize (
    const int ncid, const char* dimname, size_t* dimsize);

/* Read 3d (lev,lat,lon) slice from 4d NetCDF array (Time,lev,lat,lon)
 */
extern PetscErrorCode read3D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v);

/* Read 2d (lat,lon) slice from 3d NetCDF array (Time,lat,lon) */
extern PetscErrorCode read2D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v);

/* Read 2d (lat,lon) slice from 4d NetCDF array (Time,lat,lon) */

extern PetscErrorCode file_read_attribute (
    const int ncid, const char* name, PetscScalar* attr);
extern PetscErrorCode file_read_int_attribute (
    const int ncid, const char* name, int* attr);

extern PetscErrorCode readArray1D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    const int           n,
    PetscScalar*        a);

extern PetscErrorCode readArray (
    const int     ncid,
    const char*   varname,
    const size_t* start,
    const size_t* count,
    PetscScalar*  a);

extern PetscErrorCode write3D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v);

extern PetscErrorCode write3Ddump (const char* varname, Vec v);

#endif
