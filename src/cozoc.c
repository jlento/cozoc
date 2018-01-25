#include "context.h"
#include "defs.h"
#include "equation.h"
#include "io.h"
#include "options.h"
#include <petscdmda.h>
#include <petscksp.h>

static char help[] = BANNER
    "Solves quasi-geostrophic and generalized omega equations.\n"
    "\n"
    "Usage: mpiexec [-n procs] cozoc [-f <fname>] [-h|-Q|-G]\n"
    "               [-r <s0>,<s1>]\n"
    "\n"
    "Input file:\n"
    "  -f <fname>    Input file, NetCDF4/HDF5 format, from WRF simulation\n"
    "Which time steps to process:\n"
    "  -r <s0>,<s1>  Range of steps to compute, counting from zero\n"
    "Mode:\n"
    "  -Q            Disable quasi-geostrophic omega eq calculation\n"
    "  -G            Disable generalized omega eqs calculations\n\n";

int main (int argc, char *argv[]) {

    PetscInitialize (&argc, &argv, 0, help);

    const Options   options = new_options ();
    const NCFile    ncfile  = new_file (options);
    const Equations eqs     = new_equations (options, ncfile);
    Context         ctx     = new_context (options, ncfile);

    info (
        BANNER "Input file      : %s\n"
               "Steps in file   : %zu-%zu\n"
               "Computing steps : %zu-%zu\n\n",
        ncfile.name, 0, ctx.mt - 1, ctx.first, ctx.last);

    for (size_t istep = ctx.first; istep < ctx.last + 1; istep++) {

        info ("Step: %d\n", istep);

        update_context (ncfile, istep, ctx.first, ctx.last, &ctx);

        for (size_t ieq = 0; ieq < eqs.num_eq; ieq++) {

            Vec x = solution (eqs, ieq, ctx);

            write3D (ncfile.id, istep, eqs.id_string[ieq], x);
        }
    }

    free_context (&ctx);
    close_file (ncfile);
    PetscFinalize ();
    return 0;
}
