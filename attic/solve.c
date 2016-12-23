#include "context.h"
#include "io.h"
#include "omega.h"
#include "omegaQG.h"
#include "petscksp.h"
#include "solve.h"

#undef __FUNCT__
#define __FUNCT__ "solve"
extern PetscErrorCode solve(Context ctx) {
    KSP            ksp  = ctx->ksp;
    int            eqns = ctx->eqns;
    Vec            x;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    for (int i = 0; i < ctx->mt; i++) {
        ierr = context_update(i, ctx);
        CHKERRQ(ierr);
        for (int j = 0; j < 3; j++) {
            PetscErrorCode (*compute_operator)(KSP, Mat, Mat, void *);
            PetscErrorCode (*compute_rhs)(KSP, Vec, void *);
            char *fname;
            int   eqn = eqns & (1 << j);
            switch (eqn) {
            case OMEGA_QUASI_GEOSTROPHIC:
                compute_operator = compute_operator_omega_QG;
                compute_rhs      = compute_rhs_omega_QG;
                fname            = "ome_v_qg";
                break;
            case OMEGA_GENERALIZED:
                compute_operator = omega_compute_operator;
                compute_rhs      = omega_compute_rhs_F_V;
                fname            = "ome_v";
                break;
            default:
                continue; }

            ierr = KSPSetComputeOperators(ksp, compute_operator, ctx);
            CHKERRQ(ierr);
            ierr = KSPSetComputeRHS(ksp, compute_rhs, ctx);
            CHKERRQ(ierr);
            ierr = KSPSolve(ksp, NULL, NULL);
            CHKERRQ(ierr);
            ierr = KSPGetSolution(ksp, &x);
            CHKERRQ(ierr);
            ierr = write3D(ctx->ncid, i, fname, x);
            CHKERRQ(ierr); } }
    PetscFunctionReturn(0); }
