#include "context.h"
#include "daslice.h"
#include "derivatives.h"
#include "io.h"
#include "sigma_parameter.h"
#include "utils.h"
#include "vorticity.h"
#include "wrfnc.h"
#include <petscdmda.h>
#include <petscksp.h>


extern PetscErrorCode context_create (
    const int ncid, int skip, int* steps, int* flags, Context* p_ctx) {

    Context ctx;

    ctx = *p_ctx = (Context) malloc (sizeof (tContext) );

    /* Open Input/output file and read the dimensions */

    file_get_dimsize (ncid, dimnames[XDIM], &ctx->mx);
    file_get_dimsize (ncid, dimnames[YDIM], &ctx->my);
    file_get_dimsize (ncid, dimnames[ZDIM], &ctx->mz);
    file_get_dimsize (ncid, dimnames[TIME], &ctx->mt);

    if (skip >= (int) ctx->mt)
        ERROR (
            "Trying to skip more timesteps than is in the input file");

    if (*steps > (int) ctx->mt - skip)
        *steps = (int) ctx->mt - skip;

    /* Set up the distributed 3D array layout for 1- and 2-component
     * fields */

    DMDACreate3d (
        PETSC_COMM_WORLD,
        DM_BOUNDARY_PERIODIC,
        DM_BOUNDARY_NONE,
        DM_BOUNDARY_NONE,
        DMDA_STENCIL_BOX,
        ctx->mx,
        ctx->my,
        ctx->mz,
        PETSC_DECIDE,
        PETSC_DECIDE,
        PETSC_DECIDE,
        1,
        1,
        0,
        0,
        0,
        &ctx->da);

    DMDAGetReducedDMDA (ctx->da, 2, &ctx->da2);

    /* Set up the distributed 2D array layout (xy-plane) */

    create_subdm_plane (DMDA_Z, ctx->da, &ctx->daxy);

    /* Allocate context arrays*/

    PetscMalloc1 (ctx->mz, &ctx->Pressure);
    PetscMalloc1 (ctx->my, &ctx->Coriolis_parameter);

    DMCreateGlobalVector (ctx->da, &ctx->Sigma_parameter);
    VecDuplicate (ctx->Sigma_parameter, &ctx->Vorticity);
    VecDuplicate (ctx->Sigma_parameter, &ctx->Geopotential_height);

    DMCreateGlobalVector (ctx->da2, &ctx->Horizontal_wind);

    DMCreateGlobalVector (ctx->daxy, &ctx->Surface_pressure);
    DMCreateGlobalVector (ctx->da, &ctx->Temperature);
    VecDuplicate (ctx->Temperature, &ctx->Geopotential_height);


    /* These are read here because they do not depend on time */

    /* Pressure levels (z-coordinate) */
    {
        size_t start[1] = {0 };
        size_t count[1] = {ctx->mz };
        readArray (ncid, fieldnames[z], start, count, ctx->Pressure); }

    /* Coriollis parameter is taken to be a function of latitude, only
     */
    {
        size_t start[3] = {0, 0, 0 };
        size_t count[3] = {1, ctx->my, 1 };
        readArray (
            ncid, fieldnames[F], start, count, ctx->Coriolis_parameter); }

    /* Grid spacings */
    file_read_attribute (ncid, "DX", &ctx->hx);
    file_read_attribute (ncid, "DY", &ctx->hy);
    ctx->hz =
        ctx->Pressure[1] - ctx->Pressure[0]; /* hz is negative!!! */

    return (0); }


extern PetscErrorCode context_destroy (Context* p_ctx) {

    Context ctx = *p_ctx;

    PetscFree (ctx->Pressure);
    PetscFree (ctx->Coriolis_parameter);

    VecDestroy (&ctx->Vorticity);
    VecDestroy (&ctx->Horizontal_wind);

    VecDestroy (&ctx->Surface_pressure);
    VecDestroy (&ctx->Temperature);
    VecDestroy (&ctx->Geopotential_height);

    //        ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
    DMDestroy (&ctx->da);
    free (ctx);

    return (0); }


extern PetscErrorCode context_update (
    const int ncid, const int time, Context ctx) {

    Vec tmpvec;

    ctx->time = time;

    /* Surface pressure should be read collectively (for efficiency)? */
    read2D (ncid, time, "PSFC", ctx->Surface_pressure);
    read3D (ncid, time, "TT", ctx->Temperature);
    read3D (ncid, time, "GHT", ctx->Geopotential_height);

    /* Horizontal wind is two-component field (u,v) */
    DMGetGlobalVector (ctx->da, &tmpvec);

    for (int i = 0; i < 2; i++) {
        char* name[2] = {"UU", "VV" };
        read3D (ncid, time, name[i], tmpvec);
        VecStrideScatter (
            tmpvec, i, ctx->Horizontal_wind, INSERT_VALUES); }

    DMRestoreGlobalVector (ctx->da, &tmpvec);

    /* Calculated fields */
    sigma_parameter (ctx);
    vorticity (ctx);

    return (0); }
