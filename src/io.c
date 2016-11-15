#include "io.h"
#include <strings.h>
#include "petscdmda.h"
#include "netcdf.h"
#include "netcdf_par.h"
#include "arrays.h"
#include "utils.h"
#include "loops.h"

extern PetscErrorCode ncfile_open(const char *wrfin,int *ncid)
{
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        ierr = nc_open_par(wrfin,NC_MPIIO|NC_NETCDF4|NC_WRITE,PETSC_COMM_WORLD,
                           MPI_INFO_NULL,ncid);ERR(ierr);
        PetscFunctionReturn(0);
}

extern PetscErrorCode ncfile_get_dimsize(const int ncid,const char *dimname,
                                         size_t *dimsize)
{
        int dimid;
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        ierr = nc_inq_dimid(ncid,dimname,&dimid);ERR(ierr);
        ierr = nc_inq_dimlen(ncid,dimid,dimsize);ERR(ierr);
        PetscFunctionReturn(0);
}

extern PetscErrorCode read3D(const int ncid, const unsigned long time,  \
                             const char *varname, Vec v)
{
        DM da;
        PetscInt zs,ys,xs,zm,ym,xm;
        PetscScalar ***a;
        size_t start[4], count[4];
        int id;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);

        ierr = VecGetDM(v,&da);

        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
        start[0] = time; start[1] = zs; start[2] = ys; start[3] = xs;
        count[0] = 1;    count[1] = zm; count[2] = ym; count[3] = xm;

        ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
        ierr = nc_get_vara_double(ncid,id,start,count,
                                  &a[start[1]][start[2]][start[3]]);ERR(ierr);
        ierr = DMDAVecRestoreArray(da,v,&a);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}

extern PetscErrorCode read2D(const int ncid, const unsigned long time,  \
                             const char *varname, Vec v)
{
        DM da;
        PetscInt ys,xs,ym,xm;
        PetscScalar **a;
        size_t start[3], count[3];
        int id;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);

        ierr = VecGetDM(v,&da);

        ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
        start[0] = time; start[1] = ys; start[2] = xs;
        count[0] = 1;    count[1] = ym; count[2] = xm;

        ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
        ierr = nc_get_vara_double(ncid,id,start,count,
                                  &a[start[1]][start[2]]);ERR(ierr);
        ierr = DMDAVecRestoreArray(da,v,&a);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}

/*
  #undef __FUNCT__
  #define __FUNCT__ "readArray2D"
  extern PetscErrorCode readArray2D(const int ncid, const int time,
  const char *name,Array2D a)
  {
  size_t start[] = {time,a.ys,a.xs}, count[] = {1,a.ym,a.xm};
  int id;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = nc_inq_varid(ncid,name,&id); ERR(ierr);
  ierr = nc_get_vara_double(ncid,id,start,count,a.data);ERR(ierr);
  PetscFunctionReturn(0);
  }
*/

#undef __FUNCT__
#define __FUNCT__ "write3D"
extern PetscErrorCode write3D(const int ncid, const unsigned long time, \
                              const char *varname, Vec v)
{
        DM da;
        int id;
        PetscInt zs,ys,xs,zm,ym,xm;
        PetscScalar ***a;
        size_t start[4], count[4];
        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
        ierr = nc_var_par_access(ncid,id,NC_COLLECTIVE);ERR(ierr);
        ierr = VecGetDM(v,&da);
        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
        start[0] = time; start[1] = zs; start[2] = ys; start[3] = xs;
        count[0] = 1;    count[1] = zm; count[2] = ym; count[3] = xm;
        ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
        ierr = nc_put_vara_double(ncid,id,start,count,
                                  &a[start[1]][start[2]][start[3]]);ERR(ierr);
        ierr = DMDAVecRestoreArray(da,v,&a);CHKERRQ(ierr);
        PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "write3Ddump"
extern PetscErrorCode write3Ddump(const char *varname, Vec v)
{
        int varid,ncid,ndims=3;
        char fname[256] ="";
        char *dnames[3] = {"ZDIM","YDIM","XDIM"};
        int ds[3];
        int dimids[3];

        DM da;
        PetscInt zs,ys,xs,zm,ym,xm;
        PetscScalar ***a;
        size_t start[3], count[3];
        PetscErrorCode ierr;

        PetscFunctionBeginUser;

        ierr = VecGetDM(v,&da);
        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
        ierr = DMDAGetInfo(da,0,&ds[2],&ds[1],&ds[0],0,0,0,0,0,0,0,0,0);
        CHKERRQ(ierr);

        strcat(fname,varname); strcat(fname,".nc");
        ierr = nc_create_par(fname,NC_MPIIO|NC_NETCDF4,
                             PETSC_COMM_WORLD,MPI_INFO_NULL,&ncid); ERR(ierr);

        for(int i=0; i<ndims; i++) {
                ierr = nc_def_dim(ncid,dnames[i],ds[i],&dimids[i]); ERR(ierr);
        }

        ierr = nc_def_var(ncid,varname,NC_FLOAT,ndims,dimids,&varid);
        ierr = nc_enddef(ncid); ERR(ierr);

        start[0] = zs; start[1] = ys; start[2] = xs;
        count[0] = zm; count[1] = ym; count[2] = xm;
        ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
        ierr = nc_put_vara_double(ncid,varid,start,count,
                                  &a[zs][ys][xs]);ERR(ierr);
        nc_close(ncid);

        ierr = DMDAVecRestoreArray(da,v,&a);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}
