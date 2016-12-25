#include "constants.h"
#include "context.h"
#include "daslice.h"
#include "derivatives.h"
#include "io.h"
#include "sigma_parameter.h"
#include "utils.h"
#include "vorticity.h"
#include "wrfnc.h"
#include <math.h>
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

    DMCreateGlobalVector (ctx->daxy, &ctx->Surface_pressure);

    DMCreateGlobalVector (ctx->da, &ctx->Temperature);
    VecDuplicate (ctx->Temperature, &ctx->Sigma_parameter);
    VecDuplicate (ctx->Temperature, &ctx->Vorticity);
    VecDuplicate (ctx->Temperature, &ctx->Geopotential_height);
    VecDuplicate (ctx->Temperature, &ctx->Diabatic_heating);

    DMCreateGlobalVector (ctx->da2, &ctx->Horizontal_wind);


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

    file_read_int_attribute (ncid, "CU_PHYSICS", &ctx->cu_physics);


    return (0); }


extern PetscErrorCode context_destroy (Context* p_ctx) {

    Context ctx = *p_ctx;

    PetscFree (ctx->Pressure);
    PetscFree (ctx->Coriolis_parameter);

    VecDestroy (&ctx->Surface_pressure);

    VecDestroy (&ctx->Temperature);
    VecDestroy (&ctx->Sigma_parameter);
    VecDestroy (&ctx->Vorticity);
    VecDestroy (&ctx->Geopotential_height);
    VecDestroy (&ctx->Diabatic_heating);

    VecDestroy (&ctx->Horizontal_wind);
    VecDestroy (&ctx->Friction);


    //        ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
    DMDestroy (&ctx->da);
    free (ctx);

    return (0); }


static int horizontal_wind (Context ctx, const int ncid) {

    DM  da   = ctx->da;
    Vec V    = ctx->Horizontal_wind;
    int time = ctx->time;
    Vec tmp3d;

    DMGetGlobalVector (da, &tmp3d);

    for (int i = 0; i < 2; i++) {
        char* name[2] = {"UU", "VV" };
        read3D (ncid, time, name[i], tmp3d);
        VecStrideScatter (tmp3d, i, V, INSERT_VALUES); }

    DMRestoreGlobalVector (da, &tmp3d);
    return (0); }


static int one_over_dry_air_mass_column (
    Vec mu_inv, const int ncid, Context ctx) {

    DM  daxy = ctx->daxy;
    int time = ctx->time;
    Vec tmp2d;

    DMGetGlobalVector (daxy, &tmp2d);
    read2D (ncid, time, "MU", mu_inv);
    read2D (ncid, time, "MUB", tmp2d);
    VecAXPY (mu_inv, (PetscScalar) 1.0, tmp2d);
    VecReciprocal (mu_inv);
    DMRestoreGlobalVector (ctx->da, &tmp2d);
    return (0); }


static int diabatic_heating (Context ctx, const int ncid, Vec mvec) {

    DM            da         = ctx->da;
    DM            daxy       = ctx->daxy;
    int           cu_physics = ctx->cu_physics;
    int           time       = ctx->time;
    Vec           Q          = ctx->Diabatic_heating;
    PetscScalar*  p          = ctx->Pressure;
    const double  r          = Specific_gas_constant_of_dry_air;
    const double  cp         = Specific_heat_of_dry_air;
    Vec           tmp3d;
    PetscInt      zs, ys, xs, zm, ym, xm;
    PetscScalar **ma, ***qa;

    DMGetGlobalVector (da, &tmp3d);

    switch (cu_physics) {
    case 1:
        read3D (ncid, time, "RTHCUTEN", Q);
        read3D (ncid, time, "RTHRATEN", tmp3d);
        VecAXPY (Q, (PetscScalar) 1.0, tmp3d);
        read3D (ncid, time, "RTHBLTEN", tmp3d);
        VecAXPY (Q, (PetscScalar) 1.0, tmp3d);
        break;

    default:
        printf ("Some sensible error here"); }

    DMDAVecGetArrayRead (daxy, mvec, &ma);
    DMDAVecGetArray (da, Q, &qa);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++)
                qa[k][j][i] *= ma[j][i]; } }

    DMDAVecRestoreArray (da, Q, &qa);
    DMDAVecRestoreArrayRead (daxy, mvec, &ma);

    read3D (ncid, time, "H_DIABATIC", tmp3d);
    VecAXPY (Q, (PetscScalar) 1.0, tmp3d);

    DMDAVecGetArray (da, Q, &qa);

    for (int k = zs; k < zs + zm; k++) {
        PetscScalar alev = pow (p[k] / 10000.0, r / cp);

        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++)
                qa[k][j][i] *= alev; } }

    DMDAVecRestoreArray (da, Q, &qa);

    DMRestoreGlobalVector (ctx->da, &tmp3d);

    return (0); }


extern PetscErrorCode context_update (
    const int ncid, const int time, Context ctx) {

    DM  daxy = ctx->daxy;
    Vec mu_inv;

    ctx->time = time;

    /* Field directly from WRF output */
    read2D (ncid, time, "PSFC", ctx->Surface_pressure);

    read3D (ncid, time, "TT", ctx->Temperature);
    read3D (ncid, time, "GHT", ctx->Geopotential_height);

    horizontal_wind (ctx, ncid);

    /* Calculated fields */
    DMGetGlobalVector (daxy, &mu_inv);
    one_over_dry_air_mass_column (mu_inv, ncid, ctx);

    sigma_parameter (ctx);
    vorticity (ctx);
    diabatic_heating (ctx, ncid, mu_inv);

    DMRestoreGlobalVector (daxy, &mu_inv);

    return (0); }
