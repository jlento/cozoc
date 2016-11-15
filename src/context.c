#include "context.h"
#include "petscksp.h"
#include "petscdmda.h"
#include "netcdf.h"
#include "netcdf_par.h"
#include "io.h"
#include "daslice.h"
#include "utils.h"
#include "sigma_parameter.h"
#include "coriolis_parameter.h"
#include "vorticity.h"
#include "derivatives.h"
#include "wrfnc.h"

static PetscErrorCode output_setup(Context ctx);

#undef __FUNCT__
#define __FUNCT__ "context_create"
extern PetscErrorCode context_create(const char *fname,int *eqns,
                                     Context *p_ctx)
{
        Context         ctx;
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;

        ctx = *p_ctx = (Context)malloc(sizeof(tContext));

        ctx->eqns = *eqns;

/* Open Input/output file and read the dimensions */

        ierr = ncfile_open(fname,&ctx->ncid);CHKERRQ(ierr);
        ierr = ncfile_get_dimsize(ctx->ncid,dimnames[XDIM],&ctx->mx);
        CHKERRQ(ierr);
        ierr = ncfile_get_dimsize(ctx->ncid,dimnames[YDIM],&ctx->my);
        CHKERRQ(ierr);
        ierr = ncfile_get_dimsize(ctx->ncid,dimnames[ZDIM],&ctx->mz);
        CHKERRQ(ierr);
        ierr = ncfile_get_dimsize(ctx->ncid,dimnames[TIME],&ctx->mt);
        CHKERRQ(ierr);

/* Set up the distributed 3D array layout for scalar and 2-component fields */

        ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,
                            DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                            DMDA_STENCIL_BOX,ctx->mx,ctx->my,ctx->mz,
                            PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                            1,1,0,0,0,&ctx->da);CHKERRQ(ierr);

        ierr = DMDAGetReducedDMDA(ctx->da,2,&ctx->da2);CHKERRQ(ierr);

/* Set up the distributed 2D array layout (xy-plane) */

        ierr = create_subdm_plane(DMDA_Z,ctx->da,&ctx->daxy);

/* Set up the solver */

        ierr = KSPCreate(PETSC_COMM_WORLD,&ctx->ksp);CHKERRQ(ierr);
        ierr = KSPSetDM(ctx->ksp,ctx->da);CHKERRQ(ierr);
/*
  ierr = KSPGetPC(ctx->ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCMG); CHKERRQ(ierr);
*/
        ierr = KSPSetFromOptions(ctx->ksp);CHKERRQ(ierr);

/* Allocate context arrays*/

        ierr = PetscMalloc1(ctx->mz,&ctx->Pressure);CHKERRQ(ierr);
        ierr = PetscMalloc1(ctx->my,&ctx->Coriolis_parameter);CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(ctx->da,&ctx->Sigma_parameter);
        CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(ctx->da,&ctx->Vorticity);CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(ctx->da2,&ctx->Horizontal_wind
                );CHKERRQ(ierr);

// Old

        ierr = VecCreateSeqWithArray(
                PETSC_COMM_SELF,1,ctx->my,ctx->Coriolis_parameter,&ctx->f
                ); CHKERRQ(ierr);

        //        ierr = VecCreateSeq(PETSC_COMM_SELF,ctx->mz,&ctx->p); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(ctx->daxy,&ctx->psfc); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(ctx->da,&ctx->T); CHKERRQ(ierr);
        ierr = VecDuplicate(ctx->T,&ctx->Z); CHKERRQ(ierr);




/* These are read here because they do not depend on time */

        {
                int id;
                size_t start[1] = {0};
                size_t count[1] = {ctx->mz};
                ierr = nc_inq_varid(ctx->ncid,fieldnames[Z],&id); ERR(ierr);
                ierr = nc_get_vara_double(ctx->ncid,id,start,count,
                                          ctx->Pressure); ERR(ierr);

        }
        ierr = coriolis_parameter_read(ctx); CHKERRQ(ierr);
        ierr = nc_get_att_double(ctx->ncid,NC_GLOBAL,"DX",&ctx->hx); ERR(ierr);
        ierr = nc_get_att_double(ctx->ncid,NC_GLOBAL,"DY",&ctx->hy); ERR(ierr);
        //ierr = VecGetArray(ctx->p,&p); CHKERRQ(ierr);
        ctx->hz = ctx->Pressure[1] - ctx->Pressure[0];       /* hz is negative!!! */
        //ierr = VecRestoreArray(ctx->p,&p); CHKERRQ(ierr);

        ierr = output_setup(ctx); CHKERRQ(ierr);

        PetscFunctionReturn(0);
}


extern PetscErrorCode context_destroy(Context *p_ctx)
{
        Context         ctx = *p_ctx;
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;

        ierr = PetscFree(ctx->Coriolis_parameter);CHKERRQ(ierr);
        ierr = VecDestroy(&ctx->Vorticity);CHKERRQ(ierr);
        ierr = VecDestroy(&ctx->Horizontal_wind);CHKERRQ(ierr);


        //ierr = VecDestroy(&ctx->p);CHKERRQ(ierr);
        ierr = VecDestroy(&ctx->f);CHKERRQ(ierr);
        ierr = VecDestroy(&ctx->psfc);CHKERRQ(ierr);
        ierr = VecDestroy(&ctx->T);CHKERRQ(ierr);
        ierr = VecDestroy(&ctx->Z);CHKERRQ(ierr);

        ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
        ierr = DMDestroy(&ctx->da);CHKERRQ(ierr);
        ierr = nc_close(ctx->ncid);ERR(ierr);
        free(ctx);

        PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "context_update"
        extern PetscErrorCode context_update(const int time,Context ctx)
        {
                Vec            tmpvec;
                PetscErrorCode ierr;

                PetscFunctionBeginUser;

                ctx->time = time;
                /* Surface pressure should be read collectively (for efficiency)? */
                ierr = read2D(ctx->ncid,time,"PSFC",
                              ctx->psfc); CHKERRQ(ierr);
                ierr = read3D(ctx->ncid, time, "TT", ctx->T); CHKERRQ(ierr);
                ierr = read3D(ctx->ncid, time, "GHT", ctx->Z); CHKERRQ(ierr);

                /* Horizontal wind is two-component field (u,v) */
                ierr = DMGetGlobalVector(ctx->da,&tmpvec);CHKERRQ(ierr);
                for (int i=0; i<2; i++) {
                        char *name[2] = {"UU","VV"};
                        ierr = read3D(ctx->ncid,time,name[i],tmpvec);CHKERRQ(ierr);
                        ierr = VecStrideScatter(tmpvec,i,ctx->Horizontal_wind,
                                                INSERT_VALUES);CHKERRQ(ierr);
                }
                ierr = DMRestoreGlobalVector(ctx->da,&tmpvec);CHKERRQ(ierr);

                /* Calculated fields */
                ierr = sigma_parameter(ctx);CHKERRQ(ierr);
                ierr = vorticity(ctx);CHKERRQ(ierr);

                PetscFunctionReturn(0);
        }


static PetscErrorCode output_setup(Context ctx)
{
        const int ndims = 4;
        int ds[4];

        PetscErrorCode ierr;

        PetscFunctionBeginUser;
        for(int i = 0; i < ndims; i++) {
                ierr = nc_inq_dimid(ctx->ncid,dimnames[i],&ds[i]);ERR(ierr);
        }
        ierr = nc_redef(ctx->ncid);ERR(ierr);
        if (ctx->eqns & OMEGA_QUASI_GEOSTROPHIC) {
                int varid;
                ierr = nc_def_var(ctx->ncid,"ome_v_qg",NC_FLOAT,
                                  ndims,ds,&varid);
                if (ierr != NC_ENAMEINUSE) ERR(ierr);
        }
        if (ctx->eqns & OMEGA_GENERALIZED) {
                int varid;
                ierr = nc_def_var(ctx->ncid,"ome_v",NC_FLOAT,
                                  ndims,ds,&varid);
                if (ierr != NC_ENAMEINUSE) ERR(ierr);
        }
        ierr = nc_enddef(ctx->ncid);ERR(ierr);
        PetscFunctionReturn(0);
}
