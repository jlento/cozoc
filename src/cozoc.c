static char help[] = ""
        "Solves quasi-geostrophic and generalized omega equations, and height\n"
        "tendency equations from WRF baroclinic test case data.\n"
        "\n"
        "Usage: mpiexec [-n procs] ozoc [-f <fname>] [-h] [-t] [-Q] [-G] [-Z]\n"
        "\n"
        "Input file:\n"
        "  -f <fname>  Input file, NetCDF4/HDF5 format, from WRF simulation\n"
        "Mode:\n"
        "  -t          Run identity tests, only\n"
        "  -Q          Solve quasi-geostrophic omega equation, only\n"
        "  -G          Solve generalized omega equations, only\n"
        "  -Z          Solve height tendency equations (default), and\n"
        "              generalized omega equations if they are not already\n"
        "              in the input file\n\n";

#include <petscdmda.h>
#include <petscksp.h>
#include "context.h"
#include "idtest.h"
#include "solve.h"

static PetscErrorCode command_line_options(char *fname,PetscBool *idtestonly,
                                           int *eqns);

int main(int argc,char **argv)
{
        char            fname[PETSC_MAX_PATH_LEN] = "wrf.nc4";
        PetscBool       idtestonly = PETSC_FALSE;
        int             eqns;
        Context         ctx;
        PetscErrorCode  ierr;

        ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

        ierr = command_line_options(fname,&idtestonly,&eqns);

        ierr = context_create(fname,&eqns,&ctx);CHKERRQ(ierr);

        if (idtestonly) {
                ierr = id_test(ctx);CHKERRQ(ierr);
        } else {
                ierr = solve(ctx);CHKERRQ(ierr);
        }

        ierr = context_destroy(&ctx);
        ierr = PetscFinalize();
        return 0;
}


static PetscErrorCode command_line_options(char *fname,PetscBool *idtestonly,
                                           int *eqns)
{
        struct {char *flg; char *desc; int bf;}
        eqn[] =
                {{"-Q","Calculate quasi-geostrophic omega eq.",
                  OMEGA_QUASI_GEOSTROPHIC},
                 {"-G","Calculate generalized omega eq.",
                  OMEGA_GENERALIZED},
                 {"-Z","Calculate height-tendency eq.",
                  HEIGHT_TENDENCY}};
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;

        ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,
                                 "Options for ozoc",NULL);CHKERRQ(ierr);

        PetscOptionsBool("-t","Run identity tests, only",NULL,*idtestonly,
                         idtestonly,NULL); PetscBool  set;

        PetscOptionsString("-f","Input file, NetCDF4/HDF5 format, "
                           "from WRF simulation",NULL,fname,fname,
                           PETSC_MAX_PATH_LEN,NULL);

        *eqns = 0;
        for (int i=0; i<3; i++) {
                PetscOptionsName(eqn[i].flg,eqn[i].desc,NULL,&set);
                if (set) *eqns |= eqn[i].bf;
        }
        *eqns || (*eqns = HEIGHT_TENDENCY); // the default
        ierr = PetscOptionsEnd();CHKERRQ(ierr);
        PetscFunctionReturn(0);

}
