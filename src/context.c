#include "constants.h"
#include "context.h"
#include "daslice.h"
#include "defs.h"
#include "io.h"
#include "ops.h"
#include <math.h>
#include <petscdmda.h>
#include <petscksp.h>

Context new_context (Options const options, Files const files) {
    Context   ctx;
    const int ncid = files.ncid_in;

    ctx.ncid = ncid;

    /* Read dimensions */

    ctx.mx = files.dimsize[DIM_X];
    ctx.my = files.dimsize[DIM_Y];
    ctx.mz = files.dimsize[DIM_Z];
    ctx.mt = files.dimsize[DIM_T];

    /* Set up the distributed 3D array layout for 1- and 2-component
     * fields */

    DMDACreate3d (
        PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE,
        DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, ctx.mx, ctx.my, ctx.mz,
        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, 0, 0, 0, &ctx.da);

    DMDAGetReducedDMDA (ctx.da, 2, &ctx.da2);

    /* Set up the distributed 2D array layout (xy-plane) */

    create_subdm_plane (DMDA_Z, ctx.da, &ctx.daxy);

    KSPCreate (PETSC_COMM_WORLD, &ctx.ksp);
    KSPSetDM (ctx.ksp, ctx.da);
    KSPSetFromOptions (ctx.ksp);

    /* Allocate context arrays*/

    PetscMalloc1 (ctx.mz, &ctx.Pressure);
    PetscMalloc1 (ctx.my, &ctx.Coriolis_parameter);
    PetscMalloc1 (ctx.mt, &ctx.Time_coordinate);

    DMCreateGlobalVector (ctx.daxy, &ctx.Surface_pressure);
    VecDuplicate (ctx.Surface_pressure, &ctx.One_over_dry_air_mass_column);

    DMCreateGlobalVector (ctx.da, &ctx.Temperature);
    VecDuplicate (ctx.Temperature, &ctx.Sigma_parameter);
    VecDuplicate (ctx.Temperature, &ctx.Surface_attennuation);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity);
    VecDuplicate (ctx.Temperature, &ctx.Geopotential_height);
    VecDuplicate (ctx.Temperature, &ctx.Diabatic_heating);
    VecDuplicate (ctx.Temperature, &ctx.Diabatic_heating_attennuated);
    VecDuplicate (ctx.Temperature, &ctx.Diabatic_heating_forcing);
    VecDuplicate (ctx.Temperature, &ctx.Temperature_tendency);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency);

    for (size_t i = 0; i < NUM_GENERALIZED_OMEGA_COMPONENTS; i++)
        VecDuplicate (ctx.Temperature, &ctx.omega[i]);

    DMCreateGlobalVector (ctx.da2, &ctx.Horizontal_wind);
    VecDuplicate (ctx.Horizontal_wind, &ctx.Friction);

    /* These are read here because they are constants throughout the
     * calculation */

    /* Time coordinate */
    {
        size_t start[1] = {0 };
        size_t count[1] = {ctx.mt };
        file_read_array_double (
            ncid, "time", start, count, ctx.Time_coordinate);

        for (int i = 0; i < (int) ctx.mt; i++)
            ctx.Time_coordinate[i] *= (double) 60; }

    /* Pressure levels (z-coordinate) */
    {
        size_t start[1] = {0 };
        size_t count[1] = {ctx.mz };
        file_read_array_double (ncid, "LEV", start, count, ctx.Pressure); }

    /* Coriollis parameter is taken to be a function of latitude, only
     */
    {
        size_t start[3] = {0, 0, 0 };
        size_t count[3] = {1, ctx.my, 1 };
        file_read_array_double (
            ncid, "F", start, count, ctx.Coriolis_parameter); }

    /* Grid spacings */
    file_read_attribute (ncid, "DX", &ctx.hx);
    file_read_attribute (ncid, "DY", &ctx.hy);
    ctx.hz = ctx.Pressure[1] - ctx.Pressure[0]; /* hz is negative!!! */

    ctx.first = max_of_size_t (0, options.first);
    ctx.last  = min_of_size_t (options.last, ctx.mt - 1);

    return ctx; }

Vec new_vec (Context *ctx) {
    Vec vec;
    DMCreateGlobalVector (ctx->da, &vec);
    return vec; }

int temperature (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec T, Vec Ttend,
    Context *ctx) {

    static Vec tmpvec = 0;
    static Vec Tnext  = 0;
    static Vec Tprev  = 0;

    if (!Tnext) {    // The first step
        VecDuplicate (ctx->Temperature, &Tnext);
        VecDuplicate (ctx->Temperature, &Tprev);
        VecDuplicate (ctx->Temperature, &tmpvec); }

    if (step == first) {
        if (first == 0) {
            file_read_3d (ncid, step, "T", T);
            file_read_3d (ncid, step + 1, "T", Tnext);

            VecCopy (T, Ttend);
            VecAXPY (Ttend, -1.0, Tnext);
            VecScale (Ttend, -1.0 / (double) (t[step + 1] - t[step]) ); }
        else {
            file_read_3d (ncid, step - 1, "T", Tprev);
            file_read_3d (ncid, step, "T", T);
            file_read_3d (ncid, step + 1, "T", Tnext);

            VecCopy (Tprev, Ttend);
            VecAXPY (Ttend, -1.0, Tnext);
            VecScale (Ttend, -1.0 / (double) (t[step + 1] - t[step - 1]) ); } }
    else {
        if (step == mt - 1) {
            VecCopy (T, Tprev);
            VecCopy (Tnext, T);

            VecCopy (Tprev, Ttend);
            VecAXPY (Ttend, -1.0, T);
            VecScale (Ttend, -1.0 / (double) (t[step] - t[step - 1]) ); }
        else {
            VecCopy (T, Tprev);
            VecCopy (Tnext, T);
            file_read_3d (ncid, step + 1, "T", Tnext);

            VecCopy (Tprev, Ttend);
            VecAXPY (Ttend, -1.0, Tnext);
            VecScale (Ttend, -1.0 / (double) (t[step + 1] - t[step - 1]) ); } }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Tprev);
        VecDestroy (&Tnext); }

    return (0); }

int sigma_parameter (
    DM da, PetscInt mz, PetscScalar *p, Vec Tvec, Vec sigmavec) {
    const double   R   = Specific_gas_constant_of_dry_air;
    const double   c_p = Specific_heat_of_dry_air;
    PetscInt       zs, ys, xs, zm, ym, xm;
    PetscScalar ***T;
    PetscScalar ***sigma;

    DMDAVecGetArrayRead (da, Tvec, &T);
    DMDAVecGetArray (da, sigmavec, &sigma);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                sigma[k][j][i] =
                    log (T[k][j][i]) - (R / c_p) * log (p[k] / 100000.0); } } }

    DMDAVecRestoreArray (da, sigmavec, &sigma);
    fpder (da, mz, NULL, p, sigmavec);
    DMDAVecGetArray (da, sigmavec, &sigma);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++)
                sigma[k][j][i] *= -R * T[k][j][i] / p[k]; } }

    DMDAVecRestoreArrayRead (da, Tvec, &T);
    DMDAVecRestoreArray (da, sigmavec, &sigma);

    return (0); }

/*
int sigma_parameter (
    DM da, PetscInt mz, PetscScalar *p, Vec Tvec, Vec sigmavec) {
    const double   R   = Specific_gas_constant_of_dry_air;
    const double   c_p = Specific_heat_of_dry_air;
    PetscInt       zs, ys, xs, zm, ym, xm;
    PetscScalar ***T;
    PetscScalar ***sigma;

    // first sigma holds auxiliary variable dT/dp, temporarily
    VecCopy (Tvec, sigmavec);
    fpder (da, mz, NULL, p, sigmavec);

    // Calculating sigma --- sigma on rhs is "dT/dp"
    DMDAVecGetArrayRead (da, Tvec, &T);
    DMDAVecGetArray (da, sigmavec, &sigma);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                sigma[k][j][i] =
                    R / p[k] * (R / c_p * T[k][j][i] / p[k] - sigma[k][j][i]);

                //if (sigma[k][j][i] < sigmamin)
                //  sigma[k][j][i] = sigmamin;
            }
        }
    }

    DMDAVecRestoreArrayRead (da, Tvec, &T);
    DMDAVecRestoreArray (da, sigmavec, &sigma);

    return (0);
}
*/

/*
static int horizontal_wind_and_vorticity (const int ncid, const int step, DM da,
    DM da2, size_t my, PetscScalar hx, PetscScalar hy, Vec tmpvec, Vec V,
    Vec zeta) {

    for (int i = 0; i < 2; i++) {
        char *name[2] = {"UU", "VV"};
        file_read_3d (ncid, step, name[i], tmpvec);
        VecStrideScatter (tmpvec, i, V, INSERT_VALUES);
    }

    horizontal_rotor (da, da2, my, hx, hy, V, zeta);

    return (0);
}

static int horizontal_wind_and_vorticity_and_vorticity_tendency (int ncid,
    int step, int skip, size_t mt, double *t, DM da, DM da2, size_t my,
    PetscScalar hx, PetscScalar hy, Vec *V, Vec *Vnext, Vec *Vprev, Vec *zeta,
    Vec *zetatend, Vec *zetanext) {

    Vec tmpvec;

    // Initialize temporary vector
    tmpvec = *zetatend;

    if (step == skip) {
        if (step == 0) {
            //            tmpvec = *zetatend;
            horizontal_wind_and_vorticity (
                ncid, step, da, da2, my, hx, hy, tmpvec, *V, *zeta);
            horizontal_wind_and_vorticity (
                ncid, step + 1, da, da2, my, hx, hy, tmpvec, *Vnext, *zetanext);
            VecCopy (*zeta, *zetatend);
            VecAXPY (*zetatend, -1.0, *zetanext);
            VecScale (*zetatend, -1.0 / (double)(t[step + 1] - t[step]));
        } else {
            tmpvec = *zetatend;
            horizontal_wind_and_vorticity (
                ncid, step, da, da2, my, hx, hy, tmpvec, *V, *zeta);
            horizontal_wind_and_vorticity (
                ncid, step + 1, da, da2, my, hx, hy, tmpvec, *Vnext, *zetanext);
            horizontal_wind_and_vorticity (
                ncid, step - 1, da, da2, my, hx, hy, tmpvec, *Vprev, *zetatend);
            VecAXPY (*zetatend, -1.0, *zetanext);
            VecScale (*zetatend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    } else {
        if (step == (int)mt - 1) {
            *V = *Vnext;
            *zetatend = *zeta;
            *zeta = *zetanext;
            VecAXPY (*zetatend, -1.0, *zeta);
            VecScale (*zetatend, -1.0 / (double)(t[step] - t[step - 1]));
        } else {
            *Vprev = *V;
            *V = *Vnext;
            //            *Vnext    = tmpvec;
            //            tmpvec    = *zetatend;
            *zetatend = *zeta;
            *zeta = *zetanext;
            //            *zetanext = tmpvec;
            horizontal_wind_and_vorticity (
                ncid, step + 1, da, da2, my, hx, hy, tmpvec, *Vnext, *zetanext);
            VecAXPY (*zetatend, -1.0, *zetanext);
            VecScale (*zetatend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    }

    return (0);
}
*/

static int horizontal_wind_and_vorticity (
    const int ncid, const int step, DM da, DM da2, size_t my, PetscScalar hx,
    PetscScalar hy, Vec tmpvec, Vec V, Vec zeta) {
    for (int i = 0; i < 2; i++) {
        char name[2][2] = {"U", "V" };
        file_read_3d (ncid, step, name[i], tmpvec);
        VecStrideScatter (tmpvec, i, V, INSERT_VALUES); }

    horizontal_rotor (da, da2, my, hx, hy, V, zeta);
    return (0); }

int horizontal_wind_and_vorticity_and_vorticity_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, DM da, DM da2,
    size_t my, PetscScalar hx, PetscScalar hy, Vec V, Vec zeta, Vec zetatend,
    Context *ctx) {

    static Vec tmpvec   = 0;
    static Vec Vnext    = 0;
    static Vec Vprev    = 0;
    static Vec zetanext = 0;
    static Vec zetaprev = 0;

    info("Computing horizontal_wind_and_vorticity_and_vorticity_tendency\n");
    if (!Vnext) {    // The first step
        VecDuplicate (ctx->Horizontal_wind, &Vnext);
        VecDuplicate (ctx->Horizontal_wind, &Vprev);
        VecDuplicate (ctx->Vorticity, &zetanext);
        VecDuplicate (ctx->Vorticity, &zetaprev);
        VecDuplicate (ctx->Vorticity, &tmpvec); }

    if (step == first) {
        if (first == 0) {
            horizontal_wind_and_vorticity (
                ncid, step, da, da2, my, hx, hy, tmpvec, V, zeta);
            horizontal_wind_and_vorticity (
                ncid, step + 1, da, da2, my, hx, hy, tmpvec, Vnext, zetanext);

            VecCopy (zeta, zetatend);
            VecAXPY (zetatend, -1.0, zetanext);
            VecScale (zetatend, -1.0 / (double) (t[step + 1] - t[step]) ); }
        else {
            horizontal_wind_and_vorticity (
                ncid, step - 1, da, da2, my, hx, hy, tmpvec, Vprev, zetaprev);
            horizontal_wind_and_vorticity (
                ncid, step, da, da2, my, hx, hy, tmpvec, V, zeta);
            horizontal_wind_and_vorticity (
                ncid, step + 1, da, da2, my, hx, hy, tmpvec, Vnext, zetanext);

            VecCopy (zetaprev, zetatend);
            VecAXPY (zetatend, -1.0, zetanext);
            VecScale (zetatend, -1.0 / (double) (t[step + 1] - t[step - 1]) ); } }
    else {
        if (step == mt - 1) {
            VecCopy (V, Vprev);
            VecCopy (Vnext, V);
            VecCopy (zeta, zetaprev);
            VecCopy (zetanext, zeta);

            VecCopy (zetaprev, zetatend);
            VecAXPY (zetatend, -1.0, zeta);
            VecScale (zetatend, -1.0 / (double) (t[step] - t[step - 1]) ); }
        else {
            VecCopy (V, Vprev);
            VecCopy (Vnext, V);
            VecCopy (zeta, zetaprev);
            VecCopy (zetanext, zeta);
            horizontal_wind_and_vorticity (
                ncid, step + 1, da, da2, my, hx, hy, tmpvec, Vnext, zetanext);

            VecCopy (zetaprev, zetatend);
            VecAXPY (zetatend, -1.0, zetanext);
            VecScale (zetatend, -1.0 / (double) (t[step + 1] - t[step - 1]) ); } }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Vprev);
        VecDestroy (&Vnext);
        VecDestroy (&zetaprev);
        VecDestroy (&zetanext); }

    return (0); }

int one_over_dry_air_mass_column (
    const int ncid, const int step, Context *ctx) {

    DM  daxy = ctx->daxy;
    Vec tmp2d;

    DMGetGlobalVector (daxy, &tmp2d);
    read2D (ncid, step, "MU", ctx->One_over_dry_air_mass_column);
    read2D (ncid, step, "MUB", tmp2d);
    VecAXPY (ctx->One_over_dry_air_mass_column, (PetscScalar) 1.0, tmp2d);
    VecReciprocal (ctx->One_over_dry_air_mass_column);
    DMRestoreGlobalVector (ctx->da, &tmp2d);
    return (0); }

int diabatic_heating (Context *ctx, const int ncid, const int step) {

    DM            da   = ctx->da;
    DM            daxy = ctx->daxy;
    Vec           Q    = ctx->Diabatic_heating;
    PetscScalar * p    = ctx->Pressure;
    const double  r    = Specific_gas_constant_of_dry_air;
    const double  cp   = Specific_heat_of_dry_air;
    Vec           tmp3d;
    PetscInt      zs, ys, xs, zm, ym, xm;
    PetscScalar **ma, ***qa;

    DMGetGlobalVector (da, &tmp3d);

    file_read_3d (ncid, step, "RTHCUTEN", Q);
    file_read_3d (ncid, step, "RTHRATEN", tmp3d);
    VecAXPY (Q, (PetscScalar) 1.0, tmp3d);
    file_read_3d (ncid, step, "RTHBLTEN", tmp3d);
    VecAXPY (Q, (PetscScalar) 1.0, tmp3d);

    DMDAVecGetArrayRead (daxy, ctx->One_over_dry_air_mass_column, &ma);
    DMDAVecGetArray (da, Q, &qa);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++)
                qa[k][j][i] *= ma[j][i]; } }

    DMDAVecRestoreArray (da, Q, &qa);
    DMDAVecRestoreArrayRead (daxy, ctx->One_over_dry_air_mass_column, &ma);

    file_read_3d (ncid, step, "H_DIABATIC", tmp3d);
    VecAXPY (Q, (PetscScalar) 1.0, tmp3d);

    DMDAVecGetArray (da, Q, &qa);

    info ("ALEV: ");

    for (int k = zs; k < zs + zm; k++) {
        PetscScalar alev = pow (p[k] / 100000.0, r / cp);
        info ("%g, ", alev);

        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++)
                qa[k][j][i] *= alev; } }

    info ("\n");
    DMDAVecRestoreArray (da, Q, &qa);

    DMRestoreGlobalVector (da, &tmp3d);

    return (0); }

int friction (Context *ctx, const int ncid, const int step) {

    DM            da   = ctx->da;
    DM            da2  = ctx->da2;
    DM            daxy = ctx->daxy;
    Vec           F    = ctx->Friction;
    Vec           tmp3d;
    PetscInt      zs, ys, xs, zm, ym, xm;
    PetscScalar **ma, ****fa;

    DMGetGlobalVector (da, &tmp3d);

    file_read_3d (ncid, step, "RUCUTEN", tmp3d);
    VecStrideScatter (tmp3d, 0, F, INSERT_VALUES);
    file_read_3d (ncid, step, "RUBLTEN", tmp3d);
    VecStrideScatter (tmp3d, 0, F, ADD_VALUES);
    file_read_3d (ncid, step, "RVCUTEN", tmp3d);
    VecStrideScatter (tmp3d, 1, F, INSERT_VALUES);
    file_read_3d (ncid, step, "RVBLTEN", tmp3d);
    VecStrideScatter (tmp3d, 1, F, ADD_VALUES);

    DMDAVecGetArrayRead (daxy, ctx->One_over_dry_air_mass_column, &ma);
    DMDAVecGetArrayDOF (da2, F, &fa);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                fa[k][j][i][0] *= ma[j][i];
                fa[k][j][i][1] *= ma[j][i]; } } }

    DMDAVecRestoreArrayDOF (da2, F, &fa);
    DMDAVecRestoreArrayRead (daxy, ctx->One_over_dry_air_mass_column, &ma);

    DMRestoreGlobalVector (da, &tmp3d);

    return (0); }

/*
int context_update (const int ncid, const int step, Context *ctx) {

    DM da = ctx->da;
    DM da2 = ctx->da2;
    DM daxy = ctx->daxy;
    size_t my = ctx->my;
    size_t mz = ctx->mz;
    size_t mt = ctx->mt;
    PetscScalar hx = ctx->hx;
    PetscScalar hy = ctx->hy;
    double *time = ctx->Time_coordinate;
    int skip = ctx->skip;
    int steps = ctx->steps;
    PetscScalar *p = ctx->Pressure;
    Vec psfc = ctx->Surface_pressure;
    Vec Z = ctx->Geopotential_height;
    Vec *T = &ctx->Temperature;
    Vec *Ttend = &ctx->Temperature_tendency;
    Vec sigma = ctx->Sigma_parameter;
    Vec *V = &ctx->Horizontal_wind;
    Vec *zeta = &ctx->Vorticity;
    Vec *zetatend = &ctx->Vorticity_tendency;

    static Vec Tnext = NULL;
    static Vec Vnext = NULL;
    static Vec Vprev = NULL;
    static Vec zetanext = NULL;
    static Vec mu_inv = NULL;

    if (step == skip) { // The first step
        VecDuplicate (*T, &Tnext);
        VecDuplicate (*V, &Vnext);
        VecDuplicate (*V, &Vprev);
        VecDuplicate (*zeta, &zetanext);
        DMGetGlobalVector (daxy, &mu_inv);
    }

    read2D (ncid, step, "PSFC", psfc);
    file_read_3d (ncid, step, "GHT", Z);
    temperature (ncid, step, skip, mt, time, T, Ttend, &Tnext);
    sigma_parameter (da, mz, p, *T, sigma);
    horizontal_wind_and_vorticity_and_vorticity_tendency (ncid, step, skip, mt,
        time, da, da2, my, hx, hy, V, &Vnext, &Vprev, zeta, zetatend,
        &zetanext);
    one_over_dry_air_mass_column (mu_inv, ncid, step, ctx);
    diabatic_heating (ctx, ncid, step, mu_inv);
    friction (ctx, ncid, step, mu_inv);

    if (step == skip + steps) {
        VecDestroy (&Tnext);
        VecDestroy (&Vnext);
        VecDestroy (&zetanext);
        VecDestroy (&Tnext);
        DMRestoreGlobalVector (daxy, &mu_inv);
    }

    return (0);
}

int context_destroy (Context *ctx) {

    PetscFree (ctx->Pressure);
    PetscFree (ctx->Coriolis_parameter);

    VecDestroy (&ctx->Surface_pressure);

    VecDestroy (&ctx->Temperature);
    VecDestroy (&ctx->Sigma_parameter);
    VecDestroy (&ctx->Vorticity);
    VecDestroy (&ctx->Geopotential_height);
    VecDestroy (&ctx->Diabatic_heating);
    VecDestroy (&ctx->Temperature_tendency);
    VecDestroy (&ctx->Vorticity_tendency);

    VecDestroy (&ctx->Horizontal_wind);
    VecDestroy (&ctx->Friction);

    //        ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
    DMDestroy (&ctx->da);

    return (0);
}
*/

void update_context (size_t step, Files ncfile, Context *ctx) {

    DM           da       = ctx->da;
    DM           da2      = ctx->da2;
    size_t       my       = ctx->my;
    size_t       mz       = ctx->mz;
    size_t       mt       = ctx->mt;
    PetscScalar  hx       = ctx->hx;
    PetscScalar  hy       = ctx->hy;
    double *     time     = ctx->Time_coordinate;
    PetscScalar *p        = ctx->Pressure;
    Vec          psfc     = ctx->Surface_pressure;
    Vec          Z        = ctx->Geopotential_height;
    Vec          T        = ctx->Temperature;
    Vec          Ttend    = ctx->Temperature_tendency;
    Vec          sigma    = ctx->Sigma_parameter;
    Vec          V        = ctx->Horizontal_wind;
    Vec          zeta     = ctx->Vorticity;
    Vec          zetatend = ctx->Vorticity_tendency;

    const int ncid = ctx->ncid;

    read2D (ncid, step, "PS", psfc);
    file_read_3d (ncid, step, "GHT", Z);
    temperature (ncid, step, ctx->first, mt, time, T, Ttend, ctx);
    sigma_parameter (da, mz, p, T, sigma);
    horizontal_wind_and_vorticity_and_vorticity_tendency (
        ncid, step, ctx->first, mt, time, da, da2, my, hx, hy, V, zeta,
        zetatend, ctx);
    one_over_dry_air_mass_column (ncid, step, ctx);
    diabatic_heating (ctx, ncid, step);
    friction (ctx, ncid, step); }

void free_context (Context *ctx) {
    PetscFree (ctx->Pressure);
    PetscFree (ctx->Coriolis_parameter);

    VecDestroy (&ctx->Surface_pressure);

    VecDestroy (&ctx->Temperature);
    VecDestroy (&ctx->Sigma_parameter);
    VecDestroy (&ctx->Vorticity);
    VecDestroy (&ctx->Geopotential_height);
    VecDestroy (&ctx->Diabatic_heating);
    VecDestroy (&ctx->Temperature_tendency);
    VecDestroy (&ctx->Vorticity_tendency);

    for (size_t i = 0; i < NUM_GENERALIZED_OMEGA_COMPONENTS; i++)
        VecDestroy (&ctx->omega[i]);

    VecDestroy (&ctx->Horizontal_wind);
    VecDestroy (&ctx->Friction);

    KSPDestroy (&ctx->ksp);
    DMDestroy (&ctx->da); }
