#include "rules.h"
#include "context.h"
#include "fields.h"
#include "operators.h"
#include "targets.h"
#include <stdbool.h>

static void compute_diabatic_heating (TARGET, Targets *, Context *);

Rules new_rules (void) {
    Rules rules = {{
            [TARGET_FIELD_DIABATIC_HEATING] =
                (Rule){.prerequisites = new_target_list (TARGET_FIELD_MU_INV),
                       .recipe        = compute_diabatic_heating},

            [TARGET_FIELD_FRICTION] =
                (Rule){.prerequisites = new_target_list (TARGET_FIELD_MU_INV),
                       .recipe        = 0},

            [TARGET_FIELD_GEOPOTENTIAL_HEIGHT] =
                (Rule){.prerequisites = 0, .recipe = 0},

            [TARGET_FIELD_HORIZONTAL_WIND] =
                (Rule){.prerequisites = 0, .recipe = 0},

            [TARGET_FIELD_OMEGA_Q] = (Rule){.prerequisites = new_target_list (
                                                TARGET_FIELD_DIABATIC_HEATING,
                                                TARGET_FIELD_HORIZONTAL_WIND,
                                                TARGET_FIELD_SIGMA_PARAMETER,
                                                TARGET_FIELD_VORTICITY),
                                            .recipe = 0},

            [TARGET_FIELD_MU_INV] = (Rule){.prerequisites = 0, .recipe = 0},

            [TARGET_FIELD_TEMPERATURE] =
                (Rule){.prerequisites = 0, .recipe = 0},

            [TARGET_FIELD_TEMPERATURE_TENDENCY] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_TEMPERATURE),
                       .recipe = 0},

            [TARGET_FIELD_SIGMA_PARAMETER] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_TEMPERATURE),
                       .recipe = 0},

            [TARGET_FIELD_SURFACE_PRESSURE] =
                (Rule){.prerequisites = 0, .recipe = 0},

            [TARGET_FIELD_VORTICITY] = (Rule){.prerequisites = new_target_list (
                                                  TARGET_FIELD_HORIZONTAL_WIND),
                                              .recipe = 0},

            [TARGET_FIELD_VORTICITY_TENDENCY] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_HORIZONTAL_WIND),
                       .recipe = 0},

    }};

    return rules;
}

static char *target_name (TARGET id, const Targets *targets) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD:
        return (char *)&targets->target[id].field.name;
        break;
    case TARGET_TYPE_OPERATOR:
        return (char *)&targets->target[id].operator.name ;
        break;
    default:
        info ("ERROR: Not implemented in target_name()");
        return 0;
    }
}

void draw (const Rules *rules, const Targets *targets, const char *fname) {
    FILE *fd = fopen (fname, "w");
    PetscPrintf (
        PETSC_COMM_WORLD, "Writing target dependencies to file %s\n\n", fname);
    PetscFPrintf (PETSC_COMM_WORLD, fd, "digraph Rules {\n");
    for (size_t i = 0; i < NUM_TARGET; i++) {
        TARGETS *prereq = rules->rule[i].prerequisites;
        while (prereq) {
            PetscFPrintf (
                PETSC_COMM_WORLD, fd, "  %s -> %s\n",
                target_name (prereq->this, targets), target_name (i, targets));
            prereq = prereq->next;
        }
    }
    PetscFPrintf (PETSC_COMM_WORLD, fd, "}\n");
    fclose (fd);
}

/* Constructs a list of rules that can be executed in parallel.
 *
 * Rule is eligible for execution if
 *   (1) all prerequisites are already updated to next step, and that
 *   (2) none of the recipes of the children are pending on the current step,
 *       and
 *   (3) the rule is not at the last time step (finished).
 */

bool more_todo (
    const Rules *rules, Targets *targets, TARGETS *todo[], Context *ctx) {
    bool eligible[NUM_TARGET];
    for (size_t i = 0; i < NUM_TARGET; i++) {
        eligible[i] = true;
    }
    for (size_t i = 0; i < NUM_TARGET; i++) {
        TARGETS *prereq = rules->rule[i].prerequisites;
        while (prereq) {
            if (targets->target[i].time !=
                targets->target[prereq->this].time - 1) {    // (1)
                eligible[i] = false;
            }
            if (targets->target[i].time !=
                targets->target[prereq->this].time) {    // (2)
                eligible[prereq->this] = false;
            }
            prereq = prereq->next;
        }
        if (targets->target[i].time == ctx->last) {    // (3)
            eligible[i] = false;
        }
    }
    for (size_t i = 0; i < NUM_TARGET; i++) {
        if (eligible[i]) {
            *todo = push (i, *todo);
        }
    }
    return (*todo ? true : false);
}

void print_target_list (
    const char *title, TARGETS *head, const Targets *targets) {
    PetscPrintf (PETSC_COMM_WORLD, "%s: ", title);
    while (head) {
        info (
            "%s[%zu]%s", target_name (head->this, targets),
            targets->target[head->this].time + 1, (head->next ? ", " : ""));
        head = head->next;
    }
    info ("\n");
}

static void
compute_diabatic_heating (TARGET id, Targets *targets, Context *ctx) {
    diabatic_heating (
        ctx, targets->target[id].field.ncid, targets->target[id].time);
}

static void
read_target (TARGET id, size_t time, Targets *targets, Context *ctx) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD: {
        Field *f = &targets->target[id].field;
        file_read_3d (f->ncid, time, f->name, f->vec);
        break;
    }
    case TARGET_TYPE_OPERATOR: {
        info ("Function read() is not implemented for TARGET_TYPE_OPERATOR.\n");
        break;
    }
    default:
        info ("Internal error in rules.c, function read().\n");
    }
}

/*
void recipe (RULE id, Rules *rules) {
    Rule *field = &rules->field[id];
    field->step++;
    if (field->ncid_in)
        read_field (id, rules);
    if (field->recipe)
        field->recipe (id, rules);
}
*/

void update (TARGET id, const Rules *rules, Targets *targets, Context *ctx) {
    Target *target = &targets->target[id];
    target->time++;
    switch (target->type) {
    case TARGET_TYPE_FIELD: {
        if (rules->rule[id].recipe) {
            rules->rule[id].recipe (id, targets, ctx);
            info("Computing %s[%zu]\n", target->field.name, target->time);
        }
        break;
    }
    default:
        info ("Update of TARGET_TYPE_OPERATOR not implemented yet.\n");
    }
}

void run (const Rules *rules, Targets *targets, Context *ctx) {
    TARGETS *todo = 0;
    while ((more_todo (rules, targets, &todo, ctx))) {
        while (todo) {
            print_target_list ("todo", todo, targets);
            TARGETS *head = pop (&todo);
            update (head->this, rules, targets, ctx);
            free (head);
        }
    }
}
