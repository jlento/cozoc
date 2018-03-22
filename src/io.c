#include "defs.h"
#include "io.h"
#include "omega.h"
#include "omegaQG.h"
#include <netcdf.h>
#include <netcdf_par.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsys.h>
#include <stdbool.h>
#include <strings.h>

static bool file_exists (const char *filename) {
    FILE *fp = fopen (filename, "r");
    if (fp != 0) {
        fclose (fp);
        return true;
    } else {
        return false;
    };
}

static int io_nc_open_par (char const fname[PETSC_MAX_PATH_LEN]) {
    int ncid = -1;
    if (file_exists (fname)) {
        ERR (
            nc_open_par (
                fname, NC_MPIIO | NC_NETCDF4 | NC_WRITE, PETSC_COMM_WORLD,
                MPI_INFO_NULL, &ncid));
    } else {
        SETERRABORT (
            PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
            "Error in opening input file");
    }
    return ncid;
}

static GRIDTYPE get_grid_type (Options const *options, int const ncid) {
    return GRIDTYPE_RECTANGULAR;
}

Files new_files (const Options *options) {

    int ncidin = io_nc_open_par (options->infname);

    Files files = {
        .ncid_in   = ncidin,
        .name_in   = options->infname,
        .ncid_out  = ncidin,
        .name_out  = options->outfname,
        .grid_type = get_grid_type (options, ncidin),
        .dimname = {[DIM_T] = "time", [DIM_Z] = "vlevs",
                    [DIM_Y] = "south_north", [DIM_X] = "west_east"},
        .dimsize = {[DIM_T] = file_get_dimsize (ncidin, "time"),
                    [DIM_Z] = file_get_dimsize (ncidin, "vlevs"),
                    [DIM_Y] = file_get_dimsize (ncidin, "south_north"),
                    [DIM_X] = file_get_dimsize (ncidin, "west_east")}};

    if (options->outfname[0]) {
        int dimids[4];

        nc_create_par (
            options->outfname, NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD,
            MPI_INFO_NULL, &files.ncid_out);

        for (size_t i = 0; i < NUM_DIM; i++) {
            nc_def_dim (
                files.ncid_out, files.dimname[i], files.dimsize[i], &dimids[i]);
        }

        nc_enddef (files.ncid_out);
    }

    return files;
}

PetscErrorCode file_open (const char *wrfin, int *ncid) {
    if (file_exists (wrfin)) {
        ERR (
            nc_open_par (
                wrfin, NC_MPIIO | NC_NETCDF4 | NC_WRITE, PETSC_COMM_WORLD,
                MPI_INFO_NULL, ncid));
        return (0);
    } else {
        SETERRQ1 (
            PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Input file '%s' not found.",
            wrfin);
    }
}

void close_files (Files files) {
    nc_close (files.ncid_in);
    nc_close (files.ncid_out);
}

PetscInt file_get_dimsize (const int ncid, const char *dimname) {
    int    dimid;
    size_t dimsize;
    ERR (nc_inq_dimid (ncid, dimname, &dimid));
    ERR (nc_inq_dimlen (ncid, dimid, &dimsize));
    return dimsize;
}

PetscErrorCode
file_def_var (const int ncid, const char *name, const Files *file) {
    int dimids[NUM_DIM];

    for (int i = 0; i < NUM_DIM; i++)
        nc_inq_dimid (ncid, file->dimname[i], &dimids[i]);

    nc_def_var (ncid, name, NC_FLOAT, NUM_DIM, dimids, 0);
    return (0);
}

PetscErrorCode
file_read_attribute (const int ncid, const char *name, PetscScalar *attr) {
    nc_get_att_double (ncid, NC_GLOBAL, name, attr);
    return (0);
}

PetscErrorCode
file_read_int_attribute (const int ncid, const char *name, PetscInt *attr) {
    nc_get_att_int (ncid, NC_GLOBAL, name, attr);
    return (0);
}

int file_read_3d (
    const int ncid, const unsigned long time, const char *varname, Vec v) {

    DM             da;
    PetscScalar ***a;
    PetscInt       zs, ys, xs, zm, ym, xm;
    int            id;

    VecGetDM (v, &da);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);
    size_t start[4] = {time, zs, ys, xs};
    size_t count[4] = {1, zm, ym, xm};
    DMDAVecGetArray (da, v, &a);
    nc_inq_varid (ncid, varname, &id);
    nc_get_vara_double (
        ncid, id, start, count, &a[start[1]][start[2]][start[3]]);
    DMDAVecRestoreArray (da, v, &a);
    return (0);
}

int read2D (
    const int ncid, const unsigned long time, const char *varname, Vec v) {

    DM            da;
    PetscInt      ys, xs, ym, xm;
    PetscScalar **a;
    size_t        start[3], count[3];
    int           id;

    VecGetDM (v, &da);
    DMDAGetCorners (da, &xs, &ys, 0, &xm, &ym, 0);
    DMDAVecGetArray (da, v, &a);
    start[0] = time;
    start[1] = ys;
    start[2] = xs;
    count[0] = 1;
    count[1] = ym;
    count[2] = xm;
    nc_inq_varid (ncid, varname, &id);
    nc_get_vara_double (ncid, id, start, count, &a[start[1]][start[2]]);
    DMDAVecRestoreArray (da, v, &a);
    return (0);
}

int file_read_array_double (
    const int ncid, const char *varname, const size_t *start,
    const size_t *count, double *a_double) {
    int id;
    ERR (nc_inq_varid (ncid, varname, &id));
    nc_get_vara_double (ncid, id, start, count, a_double);
    return (0);
}

int write3D (
    const int ncid, const unsigned long time, const char *varname, Vec v) {

    DM             da;
    int            id;
    PetscInt       zs, ys, xs, zm, ym, xm;
    PetscScalar ***a;
    size_t         start[4], count[4];
    int            ierr;

    VecGetDM (v, &da);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);
    DMDAVecGetArray (da, v, &a);
    start[0] = time;
    start[1] = zs;
    start[2] = ys;
    start[3] = xs;
    count[0] = 1;
    count[1] = zm;
    count[2] = ym;
    count[3] = xm;
    info("Writing variable %s.\n", varname);
    ierr     = nc_inq_varid (ncid, varname, &id);
    ERR (ierr);
    ierr = nc_var_par_access (ncid, id, NC_COLLECTIVE);
    ERR (ierr);
    ierr = nc_put_vara_double (
        ncid, id, start, count, &a[start[1]][start[2]][start[3]]);
    ERR (ierr);
    DMDAVecRestoreArray (da, v, &a);
    return (0);
}

int write3Ddump (const char *varname, size_t mx, size_t my, size_t mz, Vec v) {
    int         varid, ncid, ndims = 4;
    char        fname[256] = "";
    const char *dnames[4]  = {"TIME", "ZDIM", "YDIM", "XDIM"};
    size_t      ds[4]      = {1, mz, my, mx};
    int         dimids[4];

    strcat (fname, varname);
    strcat (fname, ".nc");
    nc_create_par (
        fname, NC_MPIIO | NC_NETCDF4, PETSC_COMM_WORLD, MPI_INFO_NULL, &ncid);

    for (int i = 0; i < ndims; i++)
        nc_def_dim (ncid, dnames[i], ds[i], &dimids[i]);

    nc_def_var (ncid, varname, NC_FLOAT, ndims, dimids, &varid);
    nc_enddef (ncid);
    write3D (ncid, (const unsigned long)0, varname, v);
    nc_close (ncid);
    return (0);
}

/*
void write_fields (
    char const filename[], FIELD const select[], size_t const n,
    Field const fields[]) {
    if (file_exists (filename)) {

    } else {
    }
    for (size_t i = 0; i < n; i++) {
    }
}
*/
