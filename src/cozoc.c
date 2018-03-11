#include "context.h"
#include "defs.h"
#include "equation.h"
#include "io.h"
#include "options.h"
#include "rules.h"
#include "targets.h"
#include <petscsys.h>

static char help[] = BANNER
    "Solves quasi-geostrophic and generalized omega equations.\n"
    "\n"
    "Usage: [mpiexec -n procs] cozoc [-f <fname>] [-h|-Q|-G] [-r <s0>,<s1>]\n"
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

    let options = new_options ();
    let ncfile  = new_file (options);
    let rules   = new_rules ();
    var targets = new_targets (options, ncfile);

    //const Equations eqs = new_equations (options, ncfile);

    info (
        BANNER "Input file      : %s\n"
               "Steps in file   : %zu-%zu\n"
               "Computing steps : %zu-%zu\n\n",
        ncfile.name, 0, targets.context.mt - 1, targets.context.first, targets.context.last);

    draw (&rules, &targets, "deps.dot");
    run (&rules, &targets);

    /*
    Node *todo = 0;
    while (more_todo (rules, &todo)) {
        print_rule_list ("todo", todo, rules);
        Node *head = pop (&todo);
        update (head->this, &rules);
        free (head);
    }
    */

    /*
    for (size_t istep = rules.ctx.first; istep < rules.ctx.last + 1; istep++) {

        info ("Step: %d\n", istep);
        update_context (istep, ncfile, rules.ctx);

        for (size_t ieq = 0; ieq < eqs.num_eq; ieq++) {

            Vec x = solution (eqs.L[ieq], eqs.a[ieq], rules.ctx);
            write3D (ncfile.id, istep, eqs.id_string[ieq], x);
        }
    }
    */
    close_file (ncfile);
    PetscFinalize ();
    return 0;
}
