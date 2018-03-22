#include "context.h"
#include "defs.h"
#include "equation.h"
#include "io.h"
#include "options.h"
#include "rules.h"
#include "targets.h"
#include <petscsys.h>
#include <stdbool.h>

static let help = BANNER
    "Solves quasi-geostrophic and generalized omega equations.\n"
    "\n"
    "Usage: [mpiexec -n procs] cozoc [-f <fname>] [-o <fname>] [-h|-Q|-G] [-r "
    "<s0>,<s1>]\n"
    "\n"
    "  -f <fname>    Input file, NetCDF4/HDF5 format, from WRF simulation\n"
    "  -o <fname>    Output file, if separate from input file, NetCDF4/HDF5 "
    "format\n"
    "  -r <s0>,<s1>  Range of steps to compute, counting from zero\n"
    "  -Q            Disable quasi-geostrophic omega eq calculation\n"
    "  -G            Disable generalized omega eqs calculations\n\n";

int main (int argc, char *argv[]) {

    PetscInitialize (&argc, &argv, 0, help);

    let options = new_options ();
    let files = new_files (&options);
    let rules   = new_rules ();
    var ctx     = new_context (options, files);
    var targets = new_targets (options, files, &ctx);

    const Equations eqs = new_equations (options);

    info (
        BANNER "Input file            : %s\n"
               "Output file           : %s\n"
               "Steps in input file   : %zu-%zu\n"
               "Computing steps       : %zu-%zu\n\n",
        files.name_in, files.name_out, 0, ctx.mt - 1, ctx.first, ctx.last);

    draw (&rules, &targets, "deps.dot");
    run (&rules, &targets, &ctx);

    // Old main loop to be replaced by the above run()
    if (false) {
    for (size_t istep = ctx.first; istep < ctx.last + 1; istep++) {

        info ("Step: %d\n", istep);
        update_context (istep, files, &ctx);

        // for (size_t ieq = 0; ieq < eqs.num_eq; ieq++) {
        for (size_t ieq = 0; ieq < 3; ieq++) {
            Vec x = solution (eqs.L[ieq], eqs.a[ieq], ctx);
            write3D (files.ncid_out, istep, eqs.id_string[ieq], x);
        }
    }
    }
    close_files (files);
    PetscFinalize ();
    return 0;
}
