#include "petscdmda.h"

extern PetscErrorCode create_subdm_plane(DMDADirection dir,DM da,DM *dap)
{
        MPI_Comm comm,mycomm;
        DMBoundaryType bm,bn,bp,b1,b2;
        PetscInt dim,M,N,P,m,n,p,dof,s,N1,N2,n1,n2,n3;
        const PetscInt *lm,*ln,*lp,*l1,*l2,*l3;
        DMDAStencilType st;
        PetscErrorCode ierr;
        ierr = DMDAGetInfo(da,&dim,&M,&N,&P,&m,&n,&p,&dof,&s,&bm,&bn,&bp,&st);
        CHKERRQ(ierr);
        ierr = DMDAGetOwnershipRanges(da,&lm,&ln,&lp); CHKERRQ(ierr);
        switch (dir) {
        case DMDA_X :
                b1 = bn; b2 = bp;
                N1 = N; N2 = P;
                n1 = n; n2 = p; n3 = m;
                l1 = ln; l2 = lp; l3 = lm;
                break;
        case DMDA_Y :
                b1 = bm; b2 = bp;
                N1 = M; N2 = P;
                n1 = m; n2 = p; n3 = n;
                l1 = lm; l2 = lp; l3 = ln;
                break;
        case DMDA_Z :
                b1 = bm; b2 = bn;
                N1 = M; N2 = N;
                n1 = m; n2 = n; n3 = p;
                l1 = lm; l2 = ln; l3 = lp;
                break;
        }
        for (int j = 0, i = 0; i < n3; i++) {
                ierr = DMDAGetProcessorSubset(da,dir,j,&comm); CHKERRQ(ierr);
                if (comm != MPI_COMM_NULL)
                        mycomm = comm;
                j += l3[i];
        }
        ierr = DMDACreate2d(mycomm,b1,b2,st,N1,N2,n1,n2,dof,s,l1,l2,dap);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
}
