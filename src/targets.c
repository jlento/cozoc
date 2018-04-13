#include "targets.h"
#include "context.h"
#include "fields.h"
#include "netcdf.h"
#include "operators.h"
#include "options.h"
#include <stdbool.h>

Targets new_targets (Options options, Files files, Context *ctx) {
    Targets targets = {
        .target = {

                [TARGET_FIELD_DIABATIC_HEATING] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = true,
                                              .name        = "Q",
                                              .description = "Diabatic heating",
                                              .units       = "K s**-1",
                                              .vec = ctx->Diabatic_heating},
                             .time = options.first - 1},

                [TARGET_FIELD_DIABATIC_HEATING_ATTENNUATED] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){.write       = true,
                                    .name        = "Qatt",
                                    .description = "Diabatic heating with "
                                                   "surface attennuation",
                                    .units = "K s**-1",
                                    .vec   = ctx->Diabatic_heating_attennuated},
                        .time = options.first - 1},

                [TARGET_FIELD_FRICTION] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "F",
                                              .description = "Friction",
                                              .units       = 0,
                                              .vec         = ctx->Friction},
                             .time = options.first - 1},

                [TARGET_FIELD_GEOPOTENTIAL_HEIGHT] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.write       = false,
                                         .name        = "Z",
                                         .description = "Geopotential height",
                                         .units       = 0,
                                         .vec = ctx->Geopotential_height},
                             .time = options.first - 1},

                [TARGET_FIELD_HORIZONTAL_WIND] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
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
                                .write       = true,
                                .name        = "cozoc_ome_q",
                                .description = "Omega due to diabatic heating",
                                .units       = "Pa s-1",
                                .vec =
                                    ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]},
                        .time = options.first - 1},

                [TARGET_FIELD_MU_INV] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){.write = false,
                                    .name  = "mu_inv",
                                    .description =
                                        "One over dry air mass column",
                                    .units = 0,
                                    .vec   = ctx->One_over_dry_air_mass_column},
                        .time = options.first - 1},

                [TARGET_FIELD_SURFACE_ATTENNUATION] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){.write = true,
                                    .name  = "ATTENNUATION",
                                    .description =
                                        "Pressure level surface attennuation",
                                    .units = 0,
                                    .vec   = ctx->Surface_attennuation},
                        .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "TT",
                                              .description = "Temperature",
                                              .units       = "K",
                                              .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE_TENDENCY] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.write       = false,
                                         .name        = "dT",
                                         .description = "Temperature tendency",
                                         .units       = "K -s",
                                         .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_SIGMA_PARAMETER] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "sigma",
                                              .description = "Sigma parameter",
                                              .units       = "",
                                              .vec = ctx->Sigma_parameter},
                             .time = options.first - 1},

                [TARGET_FIELD_SURFACE_PRESSURE] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "PSFC",
                                              .description = "Surface pressure",
                                              .units       = "",
                                              .vec = ctx->Surface_pressure},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "zeta",
                                              .description = "Vorticity",
                                              .units       = "",
                                              .vec         = ctx->Vorticity},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.write       = false,
                                         .name        = "dzeta",
                                         .description = "Vorticity tendency",
                                         .units       = "",
                                         .vec = ctx->Vorticity_tendency},
                             .time = options.first - 1},
        }};

    nc_redef (files.ncid_out);
    for (size_t i = 0; i < NUM_TARGET; i++) {
        Target *t = &targets.target[i];
        if (t->type == TARGET_TYPE_FIELD && t->field.write) {
            file_def_var (files.ncid_out, t->field.name, &files);
        }
    }
    nc_enddef (files.ncid_out);

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
