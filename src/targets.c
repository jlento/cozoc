#include "targets.h"
#include "context.h"
#include "fields.h"
#include "operators.h"
#include "options.h"

Targets new_targets (Options options, NCFile ncfile) {
    Targets targets = {
        .context = new_context(options, ncfile),
        .target = {

                [TARGET_FIELD_DIABATIC_HEATING] =
                    (Target){.type         = TARGET_TYPE_FIELD,
                             .target.field = (Field){.ncid_in     = 0,
                                                     .ncid_out    = ncfile.id,
                                                     .name        = "diab",
                                                     .description = 0,
                                                     .units       = 0,
                                                     .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_H_DIABATIC] =
                    (Target){.type         = TARGET_TYPE_FIELD,
                             .target.field = (Field){.ncid_in  = ncfile.id,
                                                     .ncid_out = 0,
                                                     .name     = "H_DIABATIC",
                                                     .description = 0,
                                                     .units       = 0,
                                                     .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_OMEGA_Q] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .target.field =
                                 (Field){.ncid_in  = 0,
                                         .ncid_out = ncfile.id,
                                         .name     = "ome_q",
                                         .description =
                                             "omega due to diabatic heating",
                                         .units = "Pa s-1",
                                         .vec   = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_RTHBLTEN] =
                    (Target){.type         = TARGET_TYPE_FIELD,
                             .target.field = (Field){.ncid_in     = ncfile.id,
                                                     .ncid_out    = 0,
                                                     .name        = "RTHBLTEN",
                                                     .description = 0,
                                                     .units       = 0,
                                                     .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_RTHCUTEN] =
                    (Target){.type         = TARGET_TYPE_FIELD,
                             .target.field = (Field){.ncid_in     = ncfile.id,
                                                     .ncid_out    = 0,
                                                     .name        = "RTHCUTEN",
                                                     .description = 0,
                                                     .units       = 0,
                                                     .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_RTHRATEN] =
                    (Target){.type         = TARGET_TYPE_FIELD,
                             .target.field = (Field){.ncid_in     = ncfile.id,
                                                     .ncid_out    = 0,
                                                     .name        = "RTHRATEN",
                                                     .description = 0,
                                                     .units       = 0,
                                                     .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .target.field =
                                 (Field){.ncid_in     = ncfile.id,
                                         .ncid_out    = 0,
                                         .name        = "TT",
                                         .description = "Temperature",
                                         .units       = "K",
                                         .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_WIND_U] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .target.field =
                                 (Field){.ncid_in     = ncfile.id,
                                         .ncid_out    = 0,
                                         .name        = "UU",
                                         .description = "x-wind component",
                                         .units       = "m s-1",
                                         .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_WIND_V] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .target.field =
                                 (Field){.ncid_in     = ncfile.id,
                                         .ncid_out    = 0,
                                         .name        = "VV",
                                         .description = "y-wind component",
                                         .units       = "m s-1",
                                         .vec         = 0},
                             .time = options.first - 1},

                [TARGET_OPERATOR_GENERALIZED_OMEGA] =
                    (Target){.type           = TARGET_TYPE_OPERATOR,
                             .target.operator= (Operator){
                                 .name = "L",
                                 .description =
                                     "LHS operator of generalized omega eq.",
                                 .mat = 0},
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
