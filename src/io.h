#pragma once

#include "constants.h"
#include "grids.h"
#include "netcdf.h"
#include "options.h"
#include "petscdm.h"

#define NUM_DIM 4

/* NetCDF error handler */
#define ERRCODE 2
#define ERR(e)                                                                 \
    if (e) {                                                                   \
        fprintf (                                                              \
            stderr, "%s[%d]: %s: Error: %s\n", __FILE__, __LINE__, __func__,   \
            nc_strerror (e));                                                  \
        print_backtrace ();                                                    \
        exit (ERRCODE);                                                        \
    }

enum DIM { DIM_T, DIM_Z, DIM_Y, DIM_X };

typedef enum DIM         DIM;
typedef struct Dimension Dimension;

struct Dimension {
    char   name[NC_MAX_NAME + 1];
    size_t size; };

typedef struct Files Files;
struct Files {
    int         ncid_in;
    int         ncid_out;
    const char *name_in;
    const char *name_out;
    GRIDTYPE    grid_type;
    char        dimname[NUM_DIM][NC_MAX_NAME + 1];
    size_t      dimsize[NUM_DIM]; };

Files new_files (const Options *);
void  close_files (const Files);

int file_open (const char *wrfin, int *ncid);
int file_redef (const int ncid);
int file_enddef (const int ncid);
int file_def_var (const int ncid, const char *name, const Files *);

PetscInt file_get_dimsize (const int ncid, const char *dimname);

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
