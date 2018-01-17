#pragma once

#include "context.h"
#include "netcdf.h"
#include "petscdm.h"

#define DIMENSIONS TIME, ZDIM, YDIM, XDIM
#define WRF_DIMNAMES "time", "vlevs", "south_north", "west_east"

typedef enum NCFileType NCFileType;
enum NCFileType { WRF, NTYPES };

// typedef enum Dimensions Dimensions;
// enum Dimensions { DIMENSIONS, NDIMS };




enum dimensions { DIMENSIONS, NDIMS };

extern const char *dimnames[NDIMS];

enum fields { TIME_COORDINATE, Z_COORDINATE, CORIOLIS, NFIELDS };
extern const char *fieldnames[NFIELDS];

/* NetCDF error handler */
#define ERRCODE 2
#define ERR(e)                                                                 \
    if (e) {                                                                   \
        fprintf (                                                              \
            stderr, "%s[%d]: %s: Error: %s\n", __FILE__, __LINE__, __func__,   \
            nc_strerror (e));                                                  \
        exit (ERRCODE);                                                        \
    }

typedef struct NCFile NCFile;
struct NCFile {
    char name[PETSC_MAX_PATH_LEN];
    int id;
    NCFileType type;
    char dimname[NDIMS][NC_MAX_NAME + 1];
    size_t dim[NDIMS];
};

NCFile new_file (Options);

int file_open (const char *wrfin, int *ncid);
int file_close (const int ncid);
int file_redef (const int ncid);
int file_enddef (const int ncid);
int file_def_var (const int ncid, const char *name);
int file_get_dimsize (const int ncid, const char *dimname, size_t *dimsize);

/* Read 3d (lev,lat,lon) slice from 4d NetCDF array (Time,lev,lat,lon)
 */
int file_read_3d (
    const int ncid, const unsigned long time, const char *varname, Vec v);

/* Read 2d (lat,lon) slice from 3d NetCDF array (Time,lat,lon) */
int read2D (
    const int ncid, const unsigned long time, const char *varname, Vec v);

int file_read_attribute (const int ncid, const char *name, PetscScalar *attr);
int file_read_int_attribute (const int ncid, const char *name, int *attr);

int file_read_array_double (
    const int ncid, const char *varname, const size_t *start,
    const size_t *count, double *a_double);

int write3D (
    const int ncid, const unsigned long time, const char *varname, Vec v);

int write3Ddump (
    const char *varname, const size_t mx, const size_t my, const size_t mz,
    Vec v);
