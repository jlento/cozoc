#include "io.h"
#include "netcdf.h"
#include "netcdf_par.h"
#include "petscdmda.h"
#include <strings.h>


const char* dimnames[NDIMS] = {
    "time", "vlevs", "south_north", "west_east" };

const char* fieldnames[NFIELDS] = {"XTIME", "LEV", "F" };


PetscErrorCode file_open (const char* wrfin, int* ncid) {
    nc_open_par (
        wrfin,
        NC_MPIIO | NC_NETCDF4 | NC_WRITE,
        PETSC_COMM_WORLD,
        MPI_INFO_NULL,
        ncid);
    return (0); }


PetscErrorCode file_close (int ncid) {
    nc_close (ncid);
    return (0); }


PetscErrorCode file_get_dimsize (
    const int ncid, const char* dimname, size_t* dimsize) {
    int dimid;
    nc_inq_dimid (ncid, dimname, &dimid);
    nc_inq_dimlen (ncid, dimid, dimsize);
    return (0); }


PetscErrorCode file_redef (const int ncid) {
    nc_redef (ncid);
    return (0); }


PetscErrorCode file_enddef (const int ncid) {
    nc_enddef (ncid);
    return (0); }


PetscErrorCode file_def_var (const int ncid, const char* name) {
    int dimids[4];
    for (int i = 0; i < NDIMS; i++)
        nc_inq_dimid (ncid, dimnames[i], &dimids[i]);
    nc_def_var (ncid, name, NC_FLOAT, NDIMS, dimids, 0);
    return (0); }


PetscErrorCode file_read_attribute (
    const int ncid, const char* name, PetscScalar* attr) {
    nc_get_att_double (ncid, NC_GLOBAL, name, attr);
    return (0); }


PetscErrorCode file_read_int_attribute (
    const int ncid, const char* name, PetscInt* attr) {
    nc_get_att_int (ncid, NC_GLOBAL, name, attr);
    return (0); }


int file_read_3d (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v) {

    DM             da;
    PetscScalar*** a;
    PetscInt zs, ys, xs, zm, ym, xm;
    int id;

    VecGetDM (v, &da);
    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);
    size_t start[4] = {time, zs, ys, xs };
    size_t count[4] = {1, zm, ym, xm };
    DMDAVecGetArray (da, v, &a);
    nc_inq_varid (ncid, varname, &id);
    nc_get_vara_double (ncid, id, start, count, &a[start[1]][start[2]][start[3]]);
    DMDAVecRestoreArray (da, v, &a);
    return (0); }


int read2D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v) {

    DM            da;
    PetscInt      ys, xs, ym, xm;
    PetscScalar** a;
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
    nc_get_vara_double (
        ncid, id, start, count, &a[start[1]][start[2]]);
    DMDAVecRestoreArray (da, v, &a);
    return (0); }


int file_read_array_double (
    const int     ncid,
    const char*   varname,
    const size_t* start,
    const size_t* count,
    double*       a_double) {
    int id;
    nc_inq_varid (ncid, varname, &id);
    nc_get_vara_double (ncid, id, start, count, a_double);
    return (0); }


int write3D (
    const int           ncid,
    const unsigned long time,
    const char*         varname,
    Vec                 v) {

    DM             da;
    int            id;
    PetscInt       zs, ys, xs, zm, ym, xm;
    PetscScalar*** a;
    size_t         start[4], count[4];

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
    nc_inq_varid (ncid, varname, &id);
    nc_var_par_access(ncid,id,NC_COLLECTIVE);
    nc_put_vara_double (
        ncid, id, start, count, &a[start[1]][start[2]][start[3]]);
    DMDAVecRestoreArray (da, v, &a);
    return (0); }


int write3Ddump (const char* varname, size_t mx, size_t my, size_t mz, Vec v) {
    int            varid, ncid, ndims = 4;
    char           fname[256] = "";
    char*          dnames[4]  = {"TIME", "ZDIM", "YDIM", "XDIM" };
    size_t         ds[4] = {1, mz, my, mx};
    int            dimids[4];

    strcat (fname, varname);
    strcat (fname, ".nc");
    nc_create_par (
        fname,
        NC_MPIIO | NC_NETCDF4,
        PETSC_COMM_WORLD,
        MPI_INFO_NULL,
        &ncid);
    for (int i = 0; i < ndims; i++)
        nc_def_dim (ncid, dnames[i], ds[i], &dimids[i]);
    nc_def_var (ncid, varname, NC_FLOAT, ndims, dimids, &varid);
    nc_enddef (ncid);
    write3D (ncid, (const unsigned long) 0, varname, v);
    nc_close (ncid);
    return (0); }
