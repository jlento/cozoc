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
#include "equation.h"
#include "io.h"
#include "options.h"
#include <petscdmda.h>
#include <petscksp.h>

int main (int argc, char *argv[]) {
    KSP     ksp;
    Vec     x;
    Context ctx;

    PetscInitialize (&argc, &argv, 0, help);

    Options   options = new_options ();
    NCFile    ncfile  = new_file (options);
    nContext  nctx    = new_context (options, ncfile);
    Equations eqs     = new_equations (options, ncfile);

    size_t first = (options.first > 0) ? options.first : 0;
    size_t last =
        (options.last + 1 < nctx.dim[TIME]) ? options.last : nctx.dim[TIME] - 1;

    PetscPrintf (
        PETSC_COMM_WORLD, "Steps in file '%s': %zu-%zu\n"
                          "Computing steps:%*s  %zu-%zu\n",
        ncfile.name, 0, nctx.dim[TIME] - 1, strlen (ncfile.name), "", first,
        last);

    context_create (ncfile, &ctx);

    KSPCreate (PETSC_COMM_WORLD, &ksp);
    KSPSetDM (ksp, ctx.da);
    KSPSetFromOptions (ksp);

    for (size_t t = first; t < last + 1; t++) {

        PetscPrintf (PETSC_COMM_WORLD, "Step: %d\n", t);

        context_update (ncfile, t, first, last, &ctx);

        for (size_t i = 0; i < eqs.num_eq; i++) {
            KSPSetComputeOperators (ksp, eqs.L[i], &ctx);
            KSPSetComputeRHS (ksp, eqs.a[i], &ctx);
            KSPSolve (ksp, 0, 0);
            KSPGetSolution (ksp, &x);
            write3D (ncfile.id, t, eqs.id_string[i], x);
        }
    }

    KSPDestroy (&ksp);
    context_destroy (&ctx);
    file_close (ncfile.id);
    PetscFinalize ();
    return 0;
}
