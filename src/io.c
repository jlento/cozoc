#include "io.h"
#include <strings.h>
#include "petscdmda.h"
#include "netcdf.h"
#include "netcdf_par.h"
#include "arrays.h"
#include "utils.h"
#include "loops.h"
#include "config.h"

extern PetscErrorCode ncfile_open(const char *wrfin,int *ncid)
{
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
#ifdef USE_PARALLEL_NETCDF
        ierr = nc_open_par(wrfin,NC_MPIIO|NC_NETCDF4|NC_WRITE,PETSC_COMM_WORLD,
                           MPI_INFO_NULL,ncid);ERR(ierr);
#else
	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	if (rank == 0) {
		ierr = nc_open(wrfin,NC_NETCDF4|NC_WRITE,ncid);ERR(ierr);
	}
#endif
        PetscFunctionReturn(0);
}


extern PetscErrorCode ncfile_get_dimsize(const int ncid,const char *dimname,
                                         size_t *dimsize)
{
        int dimid;
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
#ifdef USE_PARALLEL_NETCDF
	ierr = nc_inq_dimid(ncid,dimname,&dimid);ERR(ierr);
	ierr = nc_inq_dimlen(ncid,dimid,dimsize);ERR(ierr);
#else
	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	if (rank == 0) {
		ierr = nc_inq_dimid(ncid,dimname,&dimid);ERR(ierr);
		ierr = nc_inq_dimlen(ncid,dimid,dimsize);ERR(ierr);
	}
	ierr = MPI_Scatter(dimsize,1,MPI_INT,dimsize,1,MPI_INT,0,
			   PETSC_COMM_WORLD);
#endif
        PetscFunctionReturn(0);
}

PetscErrorCode readAttribute(
	const int    ncid,
	const char  *name,
	PetscScalar *attr)
{
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
#ifdef USE_PARALLEL_NETCDF
        ierr = nc_get_att_double(ncid,NC_GLOBAL,name,attr); ERR(ierr);
#else
	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	if (rank == 0) {
		ierr = nc_get_att_double(ncid,NC_GLOBAL,name,attr); ERR(ierr);
	}
	ierr = MPI_Scatter(attr,1,MPI_DOUBLE,attr,1,MPI_DOUBLE,0,
			   PETSC_COMM_WORLD);
#endif
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

        ierr = VecGetDM(v,&da);

#ifndef USE_PARALLEL_NETCDF
        PetscInt mz,my,mx;
        PetscMPIInt size;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
        if (size == 1) {
#endif
                /* Parallel Netcdf4 or MPI size equal to 1 */
                ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
                start[0] = time; start[1] = zs; start[2] = ys; start[3] = xs;
                count[0] = 1;    count[1] = zm; count[2] = ym; count[3] = xm;
                ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                ierr = nc_get_vara_double(ncid,id,start,count,
                                  &a[start[1]][start[2]][start[3]]);ERR(ierr);

#ifndef USE_PARALLEL_NETCDF
        } else {
                /* Sequential Netcdf and MPI size not equal to 1 */
                PetscMPIInt rank;
                PetscScalar *buffer;
                ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

                if (rank == 0) {
                        ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);
                        CHKERRQ(ierr);
                        ierr = PetscMalloc1(mx*my*mz,&buffer);CHKERRQ(ierr);
                        start[0] = time; count[0] = 1;
                        start[1] = 0;    count[1] = mz;
                        start[2] = 0;    count[2] = my;
                        start[3] = 0;    count[3] = mx;
                        ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                        ierr = nc_get_vara_double(ncid,id,start,count,buffer);
                        ERR(ierr);

                        ierr = VecSetBlockSize(v,mx*my*mz);CHKERRQ(ierr);
                        ierr = VecSetValuesBlocked(v,1,(const PetscInt *)0,
                                           buffer,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
                        ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
                        ierr = PetscFree(buffer);CHKERRQ(ierr);
                }
        }
#endif
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

        ierr = VecGetDM(v,&da);

#ifndef USE_PARALLEL_NETCDF
        PetscInt my,mx;
        PetscMPIInt size;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
        if (size == 1) {
#endif
                /* Parallel Netcdf4 or MPI size equal to 1 */
                ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
                start[0] = time; start[1] = ys; start[2] = xs;
                count[0] = 1;    count[1] = ym; count[2] = xm;
                ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                ierr = nc_get_vara_double(ncid,id,start,count,
                                  &a[start[1]][start[2]]);ERR(ierr);

#ifndef USE_PARALLEL_NETCDF
        } else {
                /* Sequential Netcdf and MPI size not equal to 1 */
                PetscMPIInt rank;
                PetscScalar *buffer;
                ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

                if (rank == 0) {
                        ierr = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);
                        CHKERRQ(ierr);
                        ierr = PetscMalloc1(mx*my,&buffer);CHKERRQ(ierr);
                        start[0] = time; count[0] = 1;
                        start[1] = 0;    count[1] = my;
                        start[2] = 0;    count[2] = mx;
                        ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                        ierr = nc_get_vara_double(ncid,id,start,count,buffer);
                        ERR(ierr);

                        ierr = VecSetBlockSize(v,mx*my);CHKERRQ(ierr);
                        ierr = VecSetValuesBlocked(v,1,(const PetscInt *)0,
                                           buffer,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
                        ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
                        ierr = PetscFree(buffer);CHKERRQ(ierr);
                }
        }
#endif

        ierr = DMDAVecRestoreArray(da,v,&a);CHKERRQ(ierr);

        PetscFunctionReturn(0);
}


extern PetscErrorCode readArray(
        const int            ncid,
        const char          *varname,
        const size_t        *start,
        const size_t        *count,
        PetscScalar         *a)
{
        int id;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;

#ifndef USE_PARALLEL_NETCDF
        PetscMPIInt rank,size;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
        if (rank == 0) {
#endif

		/* Parallel Netcdf4 or MPI rank equal to 0 */
		ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
		ierr = nc_get_vara_double(ncid,id,start,count,a);ERR(ierr);

#ifndef USE_PARALLEL_NETCDF
	}
	/* Sequential Netcdf and MPI size > 1 */
	if (size > 1) {
		int n,dims;
		if (rank == 0) {
			ierr = nc_inq_varndims(ncid,id,&dims);ERR(ierr);
			n = count[0] - start[0];
			for (int i = 1; i < dims; i++) {
				n *= count[i] - start[i];
			}
		}
		ierr = MPI_Scatter(
			&n,1,MPI_INT,
			&n,1,MPI_INT,
			0,PETSC_COMM_WORLD);
		ierr = MPI_Scatter(
			a,n,MPI_DOUBLE,
			a,n,MPI_DOUBLE,
			0,PETSC_COMM_WORLD);
		CHKERRQ(ierr);
	}
#endif

        PetscFunctionReturn(0);
}


extern PetscErrorCode readArray1D(
        const int            ncid,
        const unsigned long  time,
        const char          *varname,
        const int            n,
        PetscScalar         *a)
{
        size_t start[2], count[2];
        int id;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;

#ifndef USE_PARALLEL_NETCDF
        PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        if (rank == 0) {
#endif
                /* Parallel Netcdf4 or MPI rank equal to 0 */
                start[0] = time; start[1] = 0;
                count[0] = 1;    count[1] = n;
                ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                ierr = nc_get_vara_double(ncid,id,start,count,a);ERR(ierr);

#ifndef USE_PARALLEL_NETCDF
        }
	ierr = MPI_Scatter(a,n,MPI_DOUBLE,a,n,MPI_DOUBLE,0,PETSC_COMM_WORLD);
        CHKERRQ(ierr);
#endif

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

        ierr = VecGetDM(v,&da);

#ifndef USE_PARALLEL_NETCDF
        PetscMPIInt size;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
        if (size == 1) {
#endif
                /* Parallel Netcdf4 or MPI size equal to 1 */
                ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,v,&a);CHKERRQ(ierr);
                start[0] = time; start[1] = zs; start[2] = ys; start[3] = xs;
                count[0] = 1;    count[1] = zm; count[2] = ym; count[3] = xm;
                ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                //ierr = nc_var_par_access(ncid,id,NC_COLLECTIVE);ERR(ierr);
                ierr = nc_put_vara_double(ncid,id,start,count,
                                          &a[start[1]][start[2]][start[3]]);
                ERR(ierr);
#ifndef USE_PARALLEL_NETCDF
        } else {
                /* Sequential Netcdf and MPI size not equal to 1 */
                PetscMPIInt rank;
                PetscInt mx,my,mz;
                Vec vloc;
                VecScatter tolocalall;
                ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
                ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);
                ierr = VecCreateSeq(PETSC_COMM_SELF,mx*my*mz,&vloc);
                CHKERRQ(ierr);
                ierr = DMDAGlobalToNaturalAllCreate(da,&tolocalall);
                CHKERRQ(ierr);
                ierr = VecScatterBegin(tolocalall,v,vloc,
                                       INSERT_VALUES,SCATTER_FORWARD);
                CHKERRQ(ierr);
                ierr = VecScatterEnd(tolocalall,v,vloc,
                                     INSERT_VALUES,SCATTER_FORWARD);
                CHKERRQ(ierr);
                ierr = DMDAVecGetArray(da,vloc,&a);CHKERRQ(ierr);
                if (rank == 0) {
                        start[0] = time; count[0] = 1;
                        start[1] = 0;    count[1] = mz;
                        start[2] = 0;    count[2] = my;
                        start[3] = 0;    count[3] = mx;
                        ierr = nc_inq_varid(ncid,varname,&id);ERR(ierr);
                        ierr = nc_put_vara_double(
                                ncid,id,start,count,
                                &a[start[1]][start[2]][start[3]]);
                        ERR(ierr);
                }
        }
#endif
        ierr = DMDAVecRestoreArray(da,v,&a);CHKERRQ(ierr);

        PetscFunctionReturn(0);

}


#undef __FUNCT__
#define __FUNCT__ "write3Ddump"
extern PetscErrorCode write3Ddump(const char *varname, Vec v)
{
        int varid,ncid,ndims=4;
        char fname[256] ="";
        char *dnames[4] = {"TIME","ZDIM","YDIM","XDIM"};
        int ds[4];
        int dimids[4];
        PetscErrorCode ierr;

        PetscFunctionBeginUser;

        strcat(fname,varname); strcat(fname,".nc");

#ifdef USE_PARALLEL_NETCDF
        ierr = nc_create_par(fname,NC_MPIIO|NC_NETCDF4,
                             PETSC_COMM_WORLD,MPI_INFO_NULL,&ncid); ERR(ierr);
#else
	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	if (rank == 0) {
                ierr = nc_create(fname,NC_NETCDF4,&ncid); ERR(ierr);
	}
#endif

#ifndef USE_PARALLEL_NETCDF
        if (rank == 0) {
#endif
                for(int i=0; i<ndims; i++) {
                        ierr = nc_def_dim(ncid,dnames[i],ds[i],&dimids[i]);
                        ERR(ierr);
                }
                ierr = nc_def_var(ncid,varname,NC_FLOAT,ndims,dimids,&varid);
                ierr = nc_enddef(ncid); ERR(ierr);
#ifndef USE_PARALLEL_NETCDF
        }
#endif

        ierr = write3D(ncid,(const unsigned long)0,varname,v);CHKERRQ(ierr);

#ifndef USE_PARALLEL_NETCDF
        if (rank == 0) {
#endif

                nc_close(ncid);

#ifndef USE_PARALLEL_NETCDF
        }
#endif

        PetscFunctionReturn(0);
}
