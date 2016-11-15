#include "coriolis_parameter.h"
#include "netcdf.h"
#include "petscdmda.h"
#include "context.h"
#include "utils.h"

/* Coriollis parameter is taken to be a function of latitude, only */

#undef __FUNCT__
#define __FUNCT__ "coriolis_parameter_read"
extern PetscErrorCode coriolis_parameter_read(Context ctx)
{
        size_t start[3] = {0, 0, 0}, count[3] = {1, ctx->my, 1};
        int id;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        ierr = nc_inq_varid(ctx->ncid, "F", &id);ERR(ierr);
        ierr = nc_get_vara_double(ctx->ncid,id,start,count,
                                  ctx->Coriolis_parameter); ERR(ierr);
        PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "coriolis_parameter_add"
extern PetscErrorCode coriolis_parameter_add(Vec x,Context ctx)
{
        PetscScalar    *f = ctx->Coriolis_parameter;
        DM             da;
        PetscInt       i,j,k,zs,ys,xs,zm,ym,xm;
        PetscScalar    ***xa;
        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        ierr = VecGetDM(x,&da);CHKERRQ(ierr);
        ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da,x,&xa);CHKERRQ(ierr);
        for (k=zs; k<zs+zm; k++) {
                for (j=ys; j<ys+ym; j++) {
                        for (i=xs; i<xs+xm; i++) {
                                xa[k][j][i] += f[j];
                        }
                }
        }
        ierr = DMDAVecRestoreArray(da,x,&xa);CHKERRQ(ierr);
        PetscFunctionReturn(0);
}
