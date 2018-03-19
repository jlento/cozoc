#include "targets.h"
#include "context.h"
#include "fields.h"
#include "operators.h"
#include "options.h"

Targets new_targets (Options options, NCFile ncfile, Context *ctx) {
    Targets targets = {
        .target = {

                [TARGET_FIELD_DIABATIC_HEATING] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = 0,
                                              .name        = "Q",
                                              .description = "Diabatic heating",
                                              .units       = 0,
                                              .vec = ctx->Diabatic_heating},
                             .time = options.first - 1},

                [TARGET_FIELD_FRICTION] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = 0,
                                              .name        = "F",
                                              .description = "Friction",
                                              .units       = 0,
                                              .vec         = ctx->Friction},
                             .time = options.first - 1},

                [TARGET_FIELD_GEOPOTENTIAL_HEIGHT] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.ncid        = 0,
                                         .name        = "Z",
                                         .description = "Geopotential height",
                                         .units       = 0,
                                         .vec = ctx->Geopotential_height},
                             .time = options.first - 1},

                [TARGET_FIELD_HORIZONTAL_WIND] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = 0,
                                              .name        = "V",
                                              .description = "Horizontal wind",
                                              .units       = 0,
                                              .vec = ctx->Horizontal_wind},
                             .time = options.first - 1},

                [TARGET_FIELD_OMEGA_Q] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){
                                .ncid        = 0,
                                .name        = "w_q",
                                .description = "Omega due to diabatic heating",
                                .units       = "Pa s-1",
                                .vec = ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]},
                        .time = options.first - 1},

                [TARGET_FIELD_MU_INV] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){.ncid = 0,
                                    .name = "mu_inv",
                                    .description =
                                        "One over dry air mass column",
                                    .units = 0,
                                    .vec   = ctx->One_over_dry_air_mass_column},
                        .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = ncfile.id,
                                              .name        = "TT",
                                              .description = "Temperature",
                                              .units       = "K",
                                              .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE_TENDENCY] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.ncid        = ncfile.id,
                                         .name        = "dT",
                                         .description = "Temperature tendency",
                                         .units       = "K -s",
                                         .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_SIGMA_PARAMETER] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = ncfile.id,
                                              .name        = "sigma",
                                              .description = "Sigma parameter",
                                              .units       = "",
                                              .vec = ctx->Sigma_parameter},
                             .time = options.first - 1},

                [TARGET_FIELD_SURFACE_PRESSURE] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = ncfile.id,
                                              .name        = "PSFC",
                                              .description = "Surface pressure",
                                              .units       = "",
                                              .vec = ctx->Surface_pressure},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.ncid        = ncfile.id,
                                              .name        = "zeta",
                                              .description = "Vorticity",
                                              .units       = "",
                                              .vec         = ctx->Vorticity},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.ncid        = ncfile.id,
                                         .name        = "dzeta",
                                         .description = "Vorticity tendency",
                                         .units       = "",
                                         .vec = ctx->Vorticity_tendency},
                             .time = options.first - 1},
        }};

    return targets;
}

TARGETS *push (TARGET target, TARGETS *oldhead) {
    TARGETS *newhead = (TARGETS *)malloc (sizeof (TARGETS));
    newhead->this    = target;
    newhead->next    = oldhead;
    return newhead;
}

TARGETS *pop (TARGETS **head) {
    TARGETS *node = *head;
    if (node) {
        *head      = node->next;
        node->next = 0;
    }
    return node;
}

TARGETS *_new_target_list (const size_t n, const TARGET f[]) {
    TARGETS *p = 0;
    for (size_t i = 0; i < n; i++) {
        p = push (f[i], p);
    }
    return p;
}
