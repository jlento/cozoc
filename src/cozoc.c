static char help[] =
    ""
    "Solves quasi-geostrophic and generalized omega equations, and height\n"
    "tendency equations from WRF baroclinic test case data.\n"
    "\n"
    "Usage: mpiexec [-n procs] ozoc [-f <fname>] [-h|-Q|-G|-Z]\n"
    "               [-s <n>] [-n <n>]\n"
    "\n"
    "Input file:\n"
    "  -f <fname>    Input file, NetCDF4/HDF5 format, from WRF simulation\n"
    "Which time steps to process:\n"
    "  -r <s1>,<s2>  Range of steps to compute, counting from zero\n"
    "Mode:\n"
    "  -Q            Disable quasi-geostrophic omega eq calculation\n"
    "  -G            Disable generalized omega eqs calculations\n\n";

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
    KSP ksp;
    Vec x;
    Context ctx;

    PetscInitialize (&argc, &argv, 0, help);

    Options options = new_options ();
    NCFile ncfile = new_file (options);

    size_t start = (options.first > 0) ? options.first : 0;
    size_t stop = (options.last + 1 < ncfile.dim[TIME]) ? (options.last + 1)
                                                        : ncfile.dim[TIME];
    PetscPrintf (
        PETSC_COMM_WORLD, "Computing %zu steps, %zu-%zu(%zu)\n",
        stop - start, start, stop - 1, ncfile.dim[TIME] - 1);

    // Context context = new_context (options, ncfile);
    PetscPrintf(PETSC_COMM_WORLD, "OK 1\n");
    context_create (ncfile.id, start, stop, &ctx);
    output_setup (ncfile.id, options);

    KSPCreate (PETSC_COMM_WORLD, &ksp);
    KSPSetDM (ksp, ctx.da);
    KSPSetFromOptions (ksp);

    for (size_t t = start; t < stop; t++) {
        PetscPrintf (PETSC_COMM_WORLD, "Time step: %d\n", t);

        context_update (ncfile.id, t, &ctx);

        if (options.compute_omega_quasi_geostrophic) {
            KSPSetComputeOperators (ksp, omega_qg_compute_operator, &ctx);
            KSPSetComputeRHS (ksp, omega_qg_compute_rhs, &ctx);
            KSPSolve (ksp, 0, 0);
            KSPGetSolution (ksp, &x);
            write3D (ncfile.id, t, OMEGA_QG_ID_STRING, x);
        }

        // write3D (ncfile.id, t, OMEGA_QG_ID_STRING, ctx->Temperature); }

        if (options.compute_omega_generalized) {
            KSPSetComputeOperators (ksp, omega_compute_operator, &ctx);

            for (int i = 0; i < N_OMEGA_COMPONENTS; i++) {
                KSPSetComputeRHS (ksp, omega_compute_rhs[i], &ctx);
                KSPSolve (ksp, 0, 0);
                KSPGetSolution (ksp, &x);
                write3D (ncfile.id, t, omega_component_id_string[i], x);
            }
        }
    }

    KSPDestroy (&ksp);
    context_destroy (&ctx);
    file_close (ncfile.id);
    PetscFinalize ();
    return 0;
}
