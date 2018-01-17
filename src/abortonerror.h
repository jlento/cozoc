#pragma once

#ifdef CHKERRQ
#undef CHKERRQ
#define CHKERRQ(n) CHKERRABORT(PETSC_COMM_SELF,(n))
#endif
