#include "rules.h"
#include "context.h"
#include "fields.h"
#include "omega.h"
#include "operators.h"
#include "targets.h"
#include <petscerror.h>
#include <petscksp.h>
#include <stdbool.h>

static void compute_diabatic_heating (TARGET, Targets *, Context *);
static void compute_friction (TARGET, Targets *, Context *);
static void compute_horizontal_wind_etc (TARGET, Targets *, Context *);
static void compute_one_over_dry_air_mass_column (TARGET, Targets *, Context *);
static void compute_temperature_and_tendency (TARGET, Targets *, Context *);
static void compute_sigma_parameter (TARGET, Targets *, Context *);
static void read_field_2d (TARGET, Targets *, Context *);
static void compute_omega_component (TARGET, Targets *, Context *);

Rules new_rules (void) {
    Rules rules = {{
            [TARGET_FIELD_DIABATIC_HEATING] =
                (Rule){.prerequisites = new_target_list (TARGET_FIELD_MU_INV),
                       .recipe        = compute_diabatic_heating},

            [TARGET_FIELD_FRICTION] =
                (Rule){.prerequisites = new_target_list (TARGET_FIELD_MU_INV),
                       .recipe        = compute_friction},

            [TARGET_FIELD_GEOPOTENTIAL_HEIGHT] =
                (Rule){.prerequisites = 0, .recipe = 0},

            [TARGET_FIELD_HORIZONTAL_WIND] =
                (Rule){.prerequisites = 0,
                       .recipe        = compute_horizontal_wind_etc},

            [TARGET_FIELD_OMEGA_Q] = (Rule){.prerequisites = new_target_list (
                                                TARGET_FIELD_DIABATIC_HEATING,
                                                TARGET_FIELD_HORIZONTAL_WIND,
                                                TARGET_FIELD_SIGMA_PARAMETER,
                                                TARGET_FIELD_VORTICITY),
                                            .recipe = compute_omega_component},

            [TARGET_FIELD_MU_INV] =
                (Rule){.prerequisites = 0,
                       .recipe        = compute_one_over_dry_air_mass_column},

            [TARGET_FIELD_TEMPERATURE] =
                (Rule){.prerequisites = 0,
                       .recipe        = compute_temperature_and_tendency},

            [TARGET_FIELD_TEMPERATURE_TENDENCY] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_TEMPERATURE),
                       .recipe = 0},

            [TARGET_FIELD_SIGMA_PARAMETER] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_TEMPERATURE),
                       .recipe = compute_sigma_parameter},

            [TARGET_FIELD_SURFACE_PRESSURE] =
                (Rule){.prerequisites = 0, .recipe = read_field_2d},

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

static void compute_friction (TARGET id, Targets *targets, Context *ctx) {
    friction (ctx, targets->target[id].field.ncid, targets->target[id].time);
}

static void
compute_omega_component (TARGET id, Targets *targets, Context *ctx) {

    KSPSetComputeOperators (ctx->ksp, omega_compute_operator, ctx);

    switch (id) {
    case TARGET_FIELD_OMEGA_Q: {
        KSPSetComputeRHS (
            ctx->ksp, omega_compute_rhs_F_Q, ctx);
        KSPSolve (ctx->ksp, 0, 0);
        KSPGetSolution (ctx->ksp, &ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]);
    }
    default: { info ("Not implemented in compute_omega_component.\n"); }
    }
}

static void
compute_horizontal_wind_etc (TARGET id, Targets *targets, Context *ctx) {
    static Vec Vnext    = NULL;
    static Vec zetanext = NULL;

    if (targets->target[id].time == ctx->first) {    // The first step
        VecDuplicate (ctx->Horizontal_wind, &Vnext);
        VecDuplicate (ctx->Vorticity, &zetanext);
    }

    horizontal_wind_and_vorticity_and_vorticity_tendency (
        targets->target[id].field.ncid, targets->target[id].time, ctx->first,
        ctx->mt, ctx->Time_coordinate, ctx->da, ctx->da2, ctx->my, ctx->hx,
        ctx->hy, &ctx->Horizontal_wind, &Vnext, &ctx->Vorticity,
        &ctx->Vorticity_tendency, &zetanext);

    if (targets->target[id].time == ctx->last) {
        PetscPrintf (
            PETSC_COMM_WORLD,
            "FIX: compute_horizontal_wind_etc fails to free vectors\n");
        /*
        VecDestroy (&Vnext);
        VecDestroy (&zetanext);
        DMRestoreGlobalVector (daxy, &mu_inv);
        */
    }
}

static void
compute_temperature_and_tendency (TARGET id, Targets *targets, Context *ctx) {
    static Vec Tnext = NULL;

    if (targets->target[id].time == ctx->first) {    // The first step
        VecDuplicate (ctx->Temperature, &Tnext);
    }

    temperature (
        targets->target[id].field.ncid, targets->target[id].time, ctx->first,
        ctx->mt, ctx->Time_coordinate, &ctx->Temperature,
        &ctx->Temperature_tendency, &Tnext);

    if (targets->target[id].time == ctx->last) {
        PetscPrintf (
            PETSC_COMM_WORLD,
            "FIX: compute_temperature fails to free vectors\n");
        /*
        VecDestroy (&Tnext);
        DMRestoreGlobalVector (daxy, &mu_inv);
        */
    }
}

static void compute_one_over_dry_air_mass_column (
    TARGET id, Targets *targets, Context *ctx) {
    one_over_dry_air_mass_column (
        targets->target[id].field.ncid, targets->target[id].time, ctx);
}

static void
compute_sigma_parameter (TARGET id, Targets *targets, Context *ctx) {
    sigma_parameter (
        ctx->da, ctx->mz, ctx->Pressure, ctx->Temperature,
        ctx->Sigma_parameter);
}

static void read_field_2d (TARGET id, Targets *targets, Context *ctx) {
    read2D (
        targets->target[id].field.ncid, targets->target[id].time,
        targets->target[id].field.name, ctx->Surface_pressure);
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
            info ("Computing %s[%zu]\n", target->field.name, target->time);
            rules->rule[id].recipe (id, targets, ctx);
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
