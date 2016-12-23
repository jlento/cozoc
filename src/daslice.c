#include <mpi.h>
#include <petscdmda.h>

PetscErrorCode create_subdm_plane (DMDADirection direction, DM da, DM *dap) {
    MPI_Comm        newcomm;
    DMBoundaryType  bm, bn, bp;
    PetscInt        dim, M, N, P, m, n, p, dof, s;
    const PetscInt *lm, *ln, *lp;
    DMDAStencilType st;
    int             color, rank;
    MPI_Comm_rank (PETSC_COMM_WORLD, &rank);
    DMDAGetInfo (da, &dim, &M, &N, &P, &m, &n, &p, &dof, &s, &bm, &bn, &bp,
                 &st);
    DMDAGetOwnershipRanges (da, &lm, &ln, &lp);

    switch (direction) {
    case DMDA_X:
        DMDAGetCorners (da, &color, 0, 0, 0, 0, 0);
        MPI_Comm_split (PETSC_COMM_WORLD, color, rank, &newcomm);
        DMDACreate2d (newcomm, bn, bp, st, N, P, n, p, dof, s, ln, lp, dap);
        break;

    case DMDA_Y:
        DMDAGetCorners (da, 0, &color, 0, 0, 0, 0);
        MPI_Comm_split (PETSC_COMM_WORLD, color, rank, &newcomm);
        DMDACreate2d (newcomm, bm, bp, st, M, P, m, p, dof, s, lm, lp, dap);
        break;

    case DMDA_Z:
        DMDAGetCorners (da, 0, 0, &color, 0, 0, 0);
        MPI_Comm_split (PETSC_COMM_WORLD, color, rank, &newcomm);
        DMDACreate2d (newcomm, bm, bn, st, M, N, m, n, dof, s, lm, ln, dap);
        break; }

    return (0); }
