static char help[] =
    ""
    "Solves quasi-geostrophic and generalized omega equations, and "
    "height\n"
    "tendency equations from WRF baroclinic test case data.\n"
    "\n"
    "Usage: mpiexec [-n procs] ozoc [-f <fname>] [-h|-Q|-G|-Z]\n"
    "               [-s <n>] [-n <n>]\n"
    "\n"
    "Input file:\n"
    "  -f <fname>  Input file, NetCDF4/HDF5 format, from WRF "
    "simulation\n"
    "Which time steps to process:\n"
    "  -s <n>      Skip <n> timesteps. Default: start from the first\n"
    "  -n <n>      Process <n> timesteps, only. Default: to the last"
    "timestep\n"
    "Mode:\n"
    "  -Q          Solve quasi-geostrophic omega equation, only\n"
    "  -G          Solve generalized omega equations, only\n\n";

    /* TODO:
    "  -Z          Solve height tendency equations (default), and\n"
    "              generalized omega equations if they are not "
    "already\n"
    "              in the input file\n\n";
    */

#include "context.h"
#include "io.h"
#include "omega.h"
#include "omegaQG.h"
#include <limits.h>
#include <petscdmda.h>
#include <petscksp.h>


static PetscErrorCode command_line_options (
    char* fname, int* skip, int* steps, int* flags) {

    PetscBool calc_Q = PETSC_FALSE;
    PetscBool calc_G = PETSC_FALSE;
    PetscBool calc_Z = PETSC_FALSE;

    PetscOptionsBegin (
        PETSC_COMM_WORLD, NULL, "Options for COZOC", NULL);

    PetscOptionsBool (
        "-Q",
        "Calculate quasi-geostrophic omega eq.",
        NULL,
        calc_Q,
        &calc_Q,
        NULL);

    PetscOptionsBool (
        "-G",
        "Calculate generalized omega eq.",
        NULL,
        calc_G,
        &calc_G,
        NULL);

    PetscOptionsBool (
        "-Z",
        "Calculate Zwack-Okossi height-tendency eq.",
        NULL,
        calc_Z,
        &calc_Z,
        NULL);

    PetscOptionsInt (
        "-s", "Skip <n> timesteps", NULL, *skip, skip, NULL);

    PetscOptionsInt (
        "-n", "Calculate <n> timesteps", NULL, *steps, steps, NULL);

    PetscOptionsString (
        "-f",
        "Input file, NetCDF4/HDF5 format, "
        "from WRF simulation",
        NULL,
        fname,
        fname,
        PETSC_MAX_PATH_LEN,
        NULL);

    *flags = 0;

    if (calc_Q)
        *flags += OMEGA_QUASI_GEOSTROPHIC;

    if (calc_G)
        *flags += OMEGA_GENERALIZED;

    if (calc_Z)
        *flags += HEIGHT_TENDENCY;

    if (*flags == 0)
        *flags = HEIGHT_TENDENCY;

    PetscOptionsEnd ();
    return (0); }


static PetscErrorCode output_setup (
    const int ncid, const int flags) {

    file_redef (ncid);

    if (flags & OMEGA_QUASI_GEOSTROPHIC) {
        file_def_var (ncid, OMEGA_QG_ID_STRING); }

    if (flags & OMEGA_GENERALIZED) {
        for (int i = 0; i < N_OMEGA_COMPONENTS; i++) {
            file_def_var (ncid, omega_component_id_string[i]); } }

    file_enddef (ncid);

    return (0); }


int main (int argc, char* argv[]) {
    char    fname[PETSC_MAX_PATH_LEN] = "wrf.nc4";
    int     flags                     = 0;
    int     skip = 0, steps = INT_MAX, ncid;
    KSP     ksp;
    Vec     x;
    Context ctx;

    PetscInitialize (&argc, &argv, NULL, help);
    command_line_options (fname, &skip, &steps, &flags);
    file_open (fname, &ncid);

    context_create (ncid, &skip, &steps, &flags, &ctx);
    output_setup (ncid, flags);

    KSPCreate (PETSC_COMM_WORLD, &ksp);
    KSPSetDM (ksp, ctx->da);
    KSPSetFromOptions (ksp);

    for (int t = skip; t < skip + steps; t++) {
        context_update (ncid, t, ctx);

        if (flags & OMEGA_QUASI_GEOSTROPHIC) {
            KSPSetComputeOperators (
                ksp, omega_qg_compute_operator, ctx);
            KSPSetComputeRHS (ksp, omega_qg_compute_rhs, ctx);
            KSPSolve (ksp, NULL, NULL);
            KSPGetSolution (ksp, &x);
            write3D (ncid, t, OMEGA_QG_ID_STRING, x); }
        //write3D (ncid, t, OMEGA_QG_ID_STRING, ctx->Temperature); }

        if (flags & OMEGA_GENERALIZED) {
            KSPSetComputeOperators (ksp, omega_compute_operator, ctx);

            for (int i = 0; i < N_OMEGA_COMPONENTS; i++) {
                KSPSetComputeRHS (ksp, omega_compute_rhs[i], ctx);
                KSPSolve (ksp, NULL, NULL);
                KSPGetSolution (ksp, &x);
                write3D (ncid, t, omega_component_id_string[i], x); } } }

    KSPDestroy (&ksp);
    context_destroy (&ctx);
    file_close (ncid);
    PetscFinalize ();
    return 0; }
