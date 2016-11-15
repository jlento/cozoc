#include "idtest.h"
#include "omegaQG.h"
#include "omega.h"
#include "petscdmda.h"
#include "io.h"
#include "utils.h"
#include "loops.h"


#undef __FUNCT__
#define __FUNCT__ "enforce_boudary_conditions"
static PetscErrorCode enforce_boudary_conditions(Vec u)
{
        DM                da;
        PetscScalar    ***a;
        PetscInt          my,mz;
        PetscErrorCode    ierr;

        PetscFunctionBeginUser;
        ierr = VecGetDM(u,&da);CHKERRQ(ierr);
        ierr = DMDAGetInfo(da,0,0,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da,u,&a);CHKERRQ(ierr);
        LOOP_KJI(da,
                 if(k == 0 || k == mz-1 || j == 0 || j == my-1)
                         a[k][j][i] = 0.0;
                );
        ierr = DMDAVecRestoreArray(da,u,&a);
        PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "create_idtest_solution_WRF_omega"
static PetscErrorCode create_idtest_solution_WRF_omega(Context ctx,Vec *u)
{
        DM da;
        PetscErrorCode       ierr;

        PetscFunctionBeginUser;
        ierr = VecGetDM(ctx->T,&da);CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(da,u);
        ierr = read3D(ctx->ncid, 0,"WW",*u);CHKERRQ(ierr);
        ierr = enforce_boudary_conditions(*u);CHKERRQ(ierr);
        PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "create_idtest_source"
static PetscErrorCode
create_idtest_source(KSP ksp, Context ctx,
                     PetscErrorCode (*compute_operator)(KSP,Mat,Mat,void*),
                     Vec u,Vec *b)
{
        DM              da;
        Mat             L;
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;
        ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
        ierr = DMCreateMatrix(da,&L);CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(da,b);CHKERRQ(ierr);
        ierr = compute_operator(ksp,NULL,L,(void*)ctx);CHKERRQ(ierr);
        ierr = MatMult(L,u,*b);CHKERRQ(ierr);
        ierr = MatDestroy(&L);CHKERRQ(ierr);
        PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "id_test"
extern PetscErrorCode id_test(Context ctx)
{
        KSP             ksp  = ctx->ksp;
        int             eqns = ctx->eqns;
        Vec             u,b,x;
        PetscErrorCode  ierr;

        PetscFunctionBeginUser;

        ierr = context_update(0,ctx);CHKERRQ(ierr);

        /* Create u,x,b */
        ierr = create_idtest_solution_WRF_omega(ctx,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(u,&x);CHKERRQ(ierr);

        for (int i = 0; i < 3; i++) {
                PetscErrorCode (*compute_operator)(KSP,Mat,Mat,void*);
                char *ename;
                int eqn = eqns & (1<<i);
                switch (eqn) {
                case OMEGA_QUASI_GEOSTROPHIC:
                        compute_operator = compute_operator_omega_QG;
                        ename = "Quasi-geostrophic omega";
                        break;
                case OMEGA_GENERALIZED:
                        compute_operator = omega_compute_operator;
                        ename = "Generelized omega";
                        break;
                default:
                        continue;
                }

                ierr = create_idtest_source(ksp,ctx,compute_operator,u,&b);
                ierr = KSPSetComputeOperators(ksp,compute_operator,
                                              ctx); CHKERRQ(ierr);

                ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

                /*
                  write3Ddump(ctx,"b",b);
                  write3Ddump(ctx,"x",x);
                */

                PetscPrintf(PETSC_COMM_WORLD,
                            "\nIdentity test\n"
                            "=============\n\n");
                PetscPrintf(PETSC_COMM_WORLD,
                            "%s equation\n"
                            "--------------------------------\n",
                            ename);

                norm2_and_absmax("   RHS source",b);
                norm2_and_absmax("Test solution",u);
                norm2_and_absmax("     Solution",x);
                ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
                norm2_and_absmax("        Error",x);
        }

        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&x); CHKERRQ(ierr);

        PetscFunctionReturn(0);
}
