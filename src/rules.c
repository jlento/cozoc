#include "context.h"
#include "equation.h"
#include "fields.h"
#include "io.h"
#include "omega.h"
#include "operators.h"
#include "ops.h"
#include "rules.h"
#include "targets.h"
#include <petscerror.h>
#include <petscksp.h>
#include <stdbool.h>

static void
compute_diabatic_heating (TARGET, Targets *, const Rules *, Context *);
static void
compute_diabatic_heating_forcing (TARGET, Targets *, const Rules *, Context *);
static void compute_friction (TARGET, Targets *, const Rules *, Context *);
static void
compute_horizontal_wind_etc (TARGET, Targets *, const Rules *, Context *);
static void
compute_omega_component (TARGET, Targets *, const Rules *, Context *);
static void compute_one_over_dry_air_mass_column (
    TARGET, Targets *, const Rules *, Context *);
static void
compute_temperature_and_tendency (TARGET, Targets *, const Rules *, Context *);
static void
compute_sigma_parameter (TARGET, Targets *, const Rules *, Context *);
static void
compute_surface_attennuation (TARGET, Targets *, const Rules *, Context *);
static void compute_surface_attennuation_factors (
    TARGET, Targets *, const Rules *, Context *);
static void read_field_2d (TARGET, Targets *, const Rules *, Context *);
static void read_field_3d (TARGET, Targets *, const Rules *, Context *);

Rules new_rules (void) {
    Rules rules = {{
            [TARGET_FIELD_DIABATIC_HEATING] =
            (Rule) {
                .prerequisites = 0,
                .recipe        = read_field_3d },

            [TARGET_FIELD_DIABATIC_HEATING_ATTENNUATED] =
            (Rule) {
                .prerequisites = new_target_list (
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_DIABATIC_HEATING),
                .recipe = compute_surface_attennuation },

            [TARGET_FIELD_DIABATIC_HEATING_FORCING] =
            (Rule) {
                .prerequisites = new_target_list (
                    TARGET_FIELD_DIABATIC_HEATING_ATTENNUATED),
                .recipe = compute_diabatic_heating_forcing },
/*
            [TARGET_FIELD_FRICTION] =
            (Rule) {
                .prerequisites = 0,
                .recipe        = compute_friction },

            [TARGET_FIELD_GEOPOTENTIAL_HEIGHT] =
            (Rule) {.prerequisites = 0, .recipe = read_field_3d },
*/

            [TARGET_FIELD_HORIZONTAL_WIND] =
            (Rule) {
                .prerequisites = 0,
                .recipe        = compute_horizontal_wind_etc },

            [TARGET_FIELD_OMEGA_Q] =
            (Rule) {
                .prerequisites = new_target_list (
                    TARGET_FIELD_DIABATIC_HEATING_FORCING,
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                .recipe = compute_omega_component },

            [TARGET_FIELD_TEMPERATURE] =
            (Rule) {
                .prerequisites = 0,
                .recipe        = compute_temperature_and_tendency },

            [TARGET_FIELD_TEMPERATURE_TENDENCY] =
            (Rule) {
                .prerequisites =
                new_target_list (TARGET_FIELD_TEMPERATURE),
                .recipe = 0 },

            [TARGET_FIELD_SIGMA_PARAMETER] =
            (Rule) {
                .prerequisites =
                new_target_list (TARGET_FIELD_TEMPERATURE),
                .recipe = compute_sigma_parameter },

            [TARGET_FIELD_SURFACE_PRESSURE] =
            (Rule) {.prerequisites = 0, .recipe = read_field_2d },

            [TARGET_FIELD_SURFACE_ATTENNUATION] =
            (Rule) {
                .prerequisites =
                new_target_list (TARGET_FIELD_SURFACE_PRESSURE),
                .recipe = compute_surface_attennuation_factors },

            [TARGET_FIELD_VORTICITY] = (Rule) {
                .prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND),
                .recipe = 0 },

            [TARGET_FIELD_VORTICITY_TENDENCY] =
            (Rule) {
                .prerequisites =
                new_target_list (TARGET_FIELD_HORIZONTAL_WIND),
                .recipe = 0 },

        } };

    return rules; }

static char *target_name (TARGET id, const Targets *targets) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD:
        return (char *) &targets->target[id].field.name;
        break;

    case TARGET_TYPE_OPERATOR:
        return (char *) &targets->target[id].operator.name ;
        break;

    default:
        info ("ERROR: Not implemented in target_name()");
        return 0; } }

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
                target_name (prereq->this, targets), target_name (i, targets) );
            prereq = prereq->next; } }

    PetscFPrintf (PETSC_COMM_WORLD, fd, "}\n");
    fclose (fd); }

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

    for (size_t i = 0; i < NUM_TARGET; i++)
        eligible[i] = true;

    for (size_t i = 0; i < NUM_TARGET; i++) {
        TARGETS *prereq = rules->rule[i].prerequisites;

        while (prereq) {
            if (targets->target[i].time !=
                    targets->target[prereq->this].time - 1)      // (1)
                eligible[i] = false;

            if (targets->target[i].time !=
                    targets->target[prereq->this].time)      // (2)
                eligible[prereq->this] = false;

            prereq = prereq->next; }

        if (targets->target[i].time == ctx->last)      // (3)
            eligible[i] = false; }

    for (size_t i = 0; i < NUM_TARGET; i++) {
        if (eligible[i])
            *todo = push (i, *todo); }

    return (*todo ? true : false); }

void print_target_list (
    const char *title, TARGETS *head, const Targets *targets) {
    PetscPrintf (PETSC_COMM_WORLD, "%s: ", title);

    while (head) {
        info (
            "%s[%zu]%s", target_name (head->this, targets),
            targets->target[head->this].time + 1, (head->next ? ", " : "") );
        head = head->next; }

    info ("\n"); }

static void compute_diabatic_heating (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    diabatic_heating (ctx, ctx->ncid, targets->target[id].time); }

static void compute_friction (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    Vec tmpvec;
    VecDuplicate (ctx->Temperature, &tmpvec);
    for (int i = 0; i < 2; i++) {
        char name[2][3] = {"FU", "FV" };
        file_read_3d (ctx->ncid, targets->target[id].time, name[i], tmpvec);
        VecStrideScatter (tmpvec, i, ctx->Friction, INSERT_VALUES); }
    VecDestroy(&tmpvec);
    }

static void compute_omega_component (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    KSPSetComputeOperators (ctx->ksp, omega_compute_operator, ctx);

    switch (id) {
    case TARGET_FIELD_OMEGA_Q: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_Q, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]);
        // KSPGetSolution (ctx->ksp, &x);
        break; }

    default: {
        info ("Not implemented in compute_omega_component.\n"); } } }

static void compute_horizontal_wind_etc (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    horizontal_wind_and_vorticity_and_vorticity_tendency (
        ctx->ncid, targets->target[id].time, ctx->first, ctx->mt,
        ctx->Time_coordinate, ctx->da, ctx->da2, ctx->my, ctx->hx, ctx->hy,
        ctx->Horizontal_wind, ctx->Vorticity, ctx->Vorticity_tendency, ctx); }

static void compute_temperature_and_tendency (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    static Vec Tnext = 0;

    if (!Tnext)      // The first step
        VecDuplicate (ctx->Temperature, &Tnext);

    temperature (
        ctx->ncid, targets->target[id].time, ctx->first, ctx->mt,
        ctx->Time_coordinate, ctx->Temperature, ctx->Temperature_tendency, ctx);

    if (targets->target[id].time == ctx->last)
        VecDestroy (&Tnext); }

static void compute_one_over_dry_air_mass_column (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    one_over_dry_air_mass_column (ctx->ncid, targets->target[id].time, ctx); }

static void compute_sigma_parameter (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    sigma_parameter (
        ctx->da, ctx->mz, ctx->Pressure, ctx->Temperature,
        ctx->Sigma_parameter); }

static void compute_surface_attennuation_factors (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    mul_fact (ctx, ctx->Surface_attennuation); }

static void compute_surface_attennuation (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    info (
        "Computing attennuation %zu, %zu\n", id,
        rules->rule[id].prerequisites->this);
    Vec y = targets->target[rules->rule[id].prerequisites->this].field.vec;
    Vec x = ctx->Surface_attennuation;
    Vec w = targets->target[id].field.vec;
    VecPointwiseMult (w, x, y); }

static void compute_diabatic_heating_forcing (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    omega_compute_rhs_F_Q (ctx->ksp, ctx->Diabatic_heating_forcing, ctx); }

static void
read_field_2d (TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    read2D (
        ctx->ncid, targets->target[id].time, targets->target[id].field.name,
        ctx->Surface_pressure); }

static void
read_field_3d (TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    file_read_3d (
        ctx->ncid, targets->target[id].time, targets->target[id].field.name,
        targets->target[id].field.vec); }

static void
read_target (TARGET id, size_t time, Targets *targets, Context *ctx) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD: {
        Field *f = &targets->target[id].field;
        file_read_3d (ctx->ncid, time, f->name, f->vec);
        break; }

    case TARGET_TYPE_OPERATOR: {
        info ("Function read() is not implemented for TARGET_TYPE_OPERATOR.\n");
        break; }

    default:
        info ("Internal error in rules.c, function read().\n"); } }

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
    Target *    target = &targets->target[id];
    const Rule *rule   = &rules->rule[id];
    target->time++;

    info ("Updating %s[%zu]\n", target->field.name, target->time);
    switch (target->type) {
    case TARGET_TYPE_FIELD: {
        if (rule->recipe) {
            rules->rule[id].recipe (id, targets, rules, ctx); }

        if (target->field.write) {
            write3D (
                ctx->ncid, target->time, target->field.name, target->field.vec); }

        break; }

    default:
        info ("Update of TARGET_TYPE_OPERATOR not implemented yet.\n"); } }

void run (const Rules *rules, Targets *targets, Context *ctx) {
    TARGETS *todo = 0;

    while ( (more_todo (rules, targets, &todo, ctx) ) ) {
        while (todo) {
            print_target_list ("todo", todo, targets);
            TARGETS *head = pop (&todo);
            update (head->this, rules, targets, ctx);
            free (head); } } }
