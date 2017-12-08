#ifndef IO_H
#define IO_H

#include "context.h"
#include "netcdf.h"
#include "petscdm.h"


enum dimensions { TIME, ZDIM, YDIM, XDIM, NDIMS };
extern const char* dimnames[NDIMS];

enum fields { TIME_COORDINATE, Z_COORDINATE, FRICTION, NFIELDS };
extern const char* fieldnames[NFIELDS];


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


int file_open (const char* wrfin, int* ncid);
int file_close (const int ncid);
int file_redef (const int ncid);
int file_enddef (const int ncid);
int file_def_var (const int ncid, const char* name);
int file_get_dimsize (
    const int ncid, const char* dimname, size_t* dimsize);

/* Read 3d (lev,lat,lon) slice from 4d NetCDF array (Time,lev,lat,lon)
 */
int file_read_3d (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v);

/* Read 2d (lat,lon) slice from 3d NetCDF array (Time,lat,lon) */
int read2D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v);

int file_read_attribute (
    const int ncid, const char* name, PetscScalar* attr);
int file_read_int_attribute (
    const int ncid, const char* name, int* attr);

int file_read_array_double (
    const int     ncid,
    const char*   varname,
    const size_t* start,
    const size_t* count,
    double*       a_double);

int write3D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v);

int write3Ddump (const char* varname, const size_t mx, const size_t my, const size_t mz, Vec v);

void handle_error(int status);

#endif
