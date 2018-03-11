#include "rules.h"
#include "context.h"
#include "fields.h"
#include "operators.h"
#include "targets.h"
#include <stdbool.h>

Rules new_rules (void) {
    Rules rules = {{

        (Rule){.target        = TARGET_FIELD_DIABATIC_HEATING,
               .prerequisites = new_target_list (
                   TARGET_FIELD_H_DIABATIC, TARGET_FIELD_RTHCUTEN,
                   TARGET_FIELD_RTHRATEN, TARGET_FIELD_RTHBLTEN),
               .recipe = 0},

        (Rule){
            .target = TARGET_FIELD_H_DIABATIC, .prerequisites = 0, .recipe = 0},

        (Rule){.target        = TARGET_FIELD_OMEGA_Q,
               .prerequisites = new_target_list (
                   TARGET_OPERATOR_GENERALIZED_OMEGA,
                   TARGET_FIELD_DIABATIC_HEATING),
               .recipe = 0},

        (Rule){
            .target = TARGET_FIELD_RTHBLTEN, .prerequisites = 0, .recipe = 0},

        (Rule){
            .target = TARGET_FIELD_RTHCUTEN, .prerequisites = 0, .recipe = 0},

        (Rule){
            .target = TARGET_FIELD_RTHRATEN, .prerequisites = 0, .recipe = 0},

        (Rule){.target        = TARGET_FIELD_TEMPERATURE,
               .prerequisites = 0,
               .recipe        = 0},

        (Rule){.target = TARGET_FIELD_WIND_U, .prerequisites = 0, .recipe = 0},

        (Rule){.target = TARGET_FIELD_WIND_V, .prerequisites = 0, .recipe = 0},

        (Rule){.target        = TARGET_OPERATOR_GENERALIZED_OMEGA,
               .prerequisites = new_target_list (
                   TARGET_FIELD_TEMPERATURE, TARGET_FIELD_WIND_U,
                   TARGET_FIELD_WIND_V),
               .recipe = 0},

    }};

    return rules;
}

static char *target_name (TARGET id, const Targets *targets) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD:
        return (char *)&targets->target[id].target.field.name;
        break;
    case TARGET_TYPE_OPERATOR:
        return (char *)&targets->target[id].target.operator.name ;
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

bool more_todo (const Rules *rules, Targets *targets, TARGETS **todo) {
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
        if (targets->target[i].time == targets->context.last) {    // (3)
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

void update (TARGET id, Targets *targets) {
    Target *target = &targets->target[id];
    target->time++;
    /*
    if (field->ncid_in)
        read_field (id, fields);
    if (field->update)
        field->update (id, fields);
    */
}

void run (const Rules *rules, Targets *targets) {
    TARGETS *todo = 0;
    while ((more_todo (rules, targets, &todo))) {
        while (todo) {
            print_target_list ("todo", todo, targets);
            TARGETS *head = pop (&todo);
            update (head->this, targets);
            free (head);
        }
    }
}
