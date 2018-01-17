static char help[] =
    ""
    "Solves quasi-geostrophic and generalized omega equations, and height\n"
    "tendency equations from WRF baroclinic test case data.\n"
    "\n"
    "Usage: mpiexec [-n procs] ozoc [-f <fname>] [-h|-Q|-G|-Z]\n"
    "               [-s <n>] [-n <n>]\n"
    "\n"
    "Input file:\n"
    "  -f <fname>  Input file, NetCDF4/HDF5 format, from WRF simulation\n"
    "Which time steps to process:\n"
    "  -s <n>      Skip <n> timesteps. Default: start from the first\n"
    "  -n <n>      Process <n> timesteps, only. Default: to the last timestep\n"
    "Mode:\n"
    "  -Q          Disable quasi-geostrophic omega eq calculation\n"
    "  -G          Disable generalized omega eqs calculations\n\n";

/* TODO:
   "  -Z          Disable height tendency eqs calculations\n\n";
*/

#include "context.h"
#include "io.h"
#include "omega.h"
#include "omegaQG.h"
#include "options.h"
#include <limits.h>
#include <petscdmda.h>
#include <petscksp.h>

static PetscErrorCode output_setup (const int ncid, const Options options) {

    file_redef (ncid);

    if (options.compute_omega_quasi_geostrophic)
        file_def_var (ncid, OMEGA_QG_ID_STRING);

    if (options.compute_omega_generalized) {
        for (int i = 0; i < N_OMEGA_COMPONENTS; i++)
            file_def_var (ncid, omega_component_id_string[i]);
    }

    file_enddef (ncid);

    return (0);
}

int main (int argc, char *argv[]) {
    int ncid;
    KSP ksp;
    Vec x;
    Context ctx;

    PetscInitialize (&argc, &argv, 0, help);

    Options options = read_options ();

    //init_context (options, &ctx);
    file_open (options.fname, &ncid);
    context_create (ncid, &options, &ctx);
    output_setup (ncid, options);



    KSPCreate (PETSC_COMM_WORLD, &ksp);
    KSPSetDM (ksp, ctx.da);
    KSPSetFromOptions (ksp);

    PetscInt skip = options.skip;
    PetscInt steps = options.steps;
    for (int t = skip; t < skip + steps; t++) {
        PetscPrintf (PETSC_COMM_WORLD, "Time step: %d\n", t);

        context_update (ncid, t, &ctx);

        if (options.compute_omega_quasi_geostrophic) {
            KSPSetComputeOperators (ksp, omega_qg_compute_operator, &ctx);
            KSPSetComputeRHS (ksp, omega_qg_compute_rhs, &ctx);
            KSPSolve (ksp, 0, 0);
            KSPGetSolution (ksp, &x);
            write3D (ncid, t, OMEGA_QG_ID_STRING, x);
        }

        // write3D (ncid, t, OMEGA_QG_ID_STRING, ctx->Temperature); }

        if (options.compute_omega_generalized) {
            KSPSetComputeOperators (ksp, omega_compute_operator, &ctx);

            for (int i = 0; i < N_OMEGA_COMPONENTS; i++) {
                KSPSetComputeRHS (ksp, omega_compute_rhs[i], &ctx);
                KSPSolve (ksp, 0, 0);
                KSPGetSolution (ksp, &x);
                write3D (ncid, t, omega_component_id_string[i], x);
            }
        }
    }

    KSPDestroy (&ksp);
    context_destroy (&ctx);
    file_close (ncid);
    PetscFinalize ();
    return 0;
}
