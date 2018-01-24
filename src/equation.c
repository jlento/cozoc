#include "equation.h"
#include "io.h"
#include "options.h"
#include "grids.h"
#include "omegaQG.h"
#include "omega.h"

Equations new_equations (Options const options, NCFile const ncfile) {
    Equations eqs = {.num_eq = 0, .L = {0}, .a = {0}, .id_string = {""}};
    if (options.compute_omega_quasi_geostrophic) {
        eqs.L[eqs.num_eq] = omega_qg_compute_operator;
        eqs.a[eqs.num_eq] = omega_qg_compute_rhs;
        strcpy(eqs.id_string[eqs.num_eq], OMEGA_QG_ID_STRING);
        eqs.num_eq++;
    }
    if (options.compute_omega_quasi_geostrophic) {
        for (size_t i = 0; i < NUM_GENERALIZED_OMEGA_COMPONENTS; i++) {
            eqs.L[eqs.num_eq] = omega_compute_operator;
            eqs.a[eqs.num_eq] = omega_compute_rhs[i];
            strcpy(eqs.id_string[eqs.num_eq],omega_component_id_string[i]);
            eqs.num_eq++;
        }
    }
    return eqs;
}
