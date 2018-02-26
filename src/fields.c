#include "fields.h"
#include "context.h"
#include "defs.h"
#include <petscsys.h>
#include <stdio.h>

static void read_field (FIELD, Fields *);

Fields new_fields (Options options, NCFile ncfile) {
    Fields fields =    //
        {.ctx   = new_context (options, ncfile),
         .field = {[FIELD_OMEGA_Q] =    //
                   (Field){.step            = options.first - 1,
                           .ncid_in         = 0,
                           .ncid_out        = ncfile.id,
                           .name            = "ome_q",
                           .description     = "omega due to diabatic heating",
                           .units           = "Pa s-1",
                           .parents         = new_list (FIELD_DIABATIC_HEATING),
                           .children        = 0,
                           .update          = 0},
                   [FIELD_DIABATIC_HEATING] =    //
                   (Field){.step        = options.first - 1,
                           .ncid_in     = 0,
                           .ncid_out    = ncfile.id,
                           .name        = "diab",
                           .description = 0,
                           .units       = 0,
                           .parents     = new_list (
                               FIELD_RTHCUTEN, FIELD_RTHRATEN, FIELD_RTHBLTEN,
                               FIELD_H_DIABATIC),
                           .children  = 0,
                           .update    = 0},
                   [FIELD_H_DIABATIC] =    //
                   (Field){.step        = options.first - 1,
                           .ncid_in     = ncfile.id,
                           .ncid_out    = 0,
                           .name        = "H_DIABATIC",
                           .description = 0,
                           .units       = 0,
                           .vec         = new_vec (&fields.ctx),
                           .parents     = 0,
                           .children    = 0,
                           .update      = 0},
                   [FIELD_RTHCUTEN]     =    //
                   (Field){.step        = options.first - 1,
                           .ncid_in     = ncfile.id,
                           .ncid_out    = 0,
                           .name        = "RTHCUTEN",
                           .description = 0,
                           .units       = 0,
                           .vec         = new_vec (&fields.ctx),
                           .parents     = 0,
                           .children    = 0,
                           .update      = 0},
                   [FIELD_RTHRATEN]     =    //
                   (Field){.step        = options.first - 1,
                           .ncid_in     = ncfile.id,
                           .ncid_out    = 0,
                           .name        = "RTHRATEN",
                           .description = 0,
                           .units       = 0,
                           .vec         = new_vec (&fields.ctx),
                           .parents     = 0,
                           .children    = 0,
                           .update      = 0},
                   [FIELD_RTHBLTEN]     =    //
                   (Field){.step        = options.first - 1,
                           .ncid_in     = ncfile.id,
                           .ncid_out    = 0,
                           .name        = "RTHBLTEN",
                           .description = 0,
                           .units       = 0,
                           .vec         = new_vec (&fields.ctx),
                           .parents     = 0,
                           .children    = 0,
                           .update      = 0}}};

    Node *parent = fields.field[0].parents;
    for (size_t i = 0; i < NUM_FIELD; i++, parent = fields.field[i].parents) {
        while (parent) {
            fields.field[parent->this].children =
                push (i, fields.field[parent->this].children);
            parent = parent->next;
        }
    }

    return fields;
}

void read_field (FIELD id, Fields *fields) {
    Field field = fields->field[id];
    file_read_3d (field.ncid_in, field.step, field.name, field.vec);
}

Node *push (FIELD field, Node *oldhead) {
    Node *newhead = (Node *)malloc (sizeof (Node));
    newhead->this = field;
    newhead->next = oldhead;
    return newhead;
}

Node *pop (Node **head) {
    Node *node = *head;
    if (node) {
        *head      = node->next;
        node->next = 0;
    }
    return node;
}

void draw_tree (Fields fields, const char *fname) {
    FILE *fd = fopen (fname, "w");
    PetscPrintf (
        PETSC_COMM_WORLD, "Writing field dependencies to file %s\n\n", fname);
    PetscFPrintf (PETSC_COMM_WORLD, fd, "digraph Fields {\n");
    Node *parent = fields.field[0].parents;
    for (size_t i = 0; i < NUM_FIELD; i++, parent = fields.field[i].parents) {
        while (parent) {
            PetscFPrintf (
                PETSC_COMM_WORLD, fd, "  %s -> %s\n",
                fields.field[parent->this].name, fields.field[i].name);
            parent = parent->next;
        }
    }
    PetscFPrintf (PETSC_COMM_WORLD, fd, "}\n");
    fclose (fd);
}

/* Constructs a list of fields that can be updated.
 *
 * Requires that
 *   (1) all parents are already updated to next step, and that
 *   (2) none of the updates of the children are pending on the current step
 */

Node *more_todo (Fields fields, Node **todo) {
    Node *parent = fields.field[0].parents;
    Node *child  = fields.field[0].children;
    for (size_t i = 0; i < NUM_FIELD; i++, parent = fields.field[i].parents,
                child = fields.field[i].children) {
        Node *p = *todo;
        if (fields.field[i].step == fields.ctx.last) {
            goto BREAK;
        }
        while (parent) {
            if (fields.field[parent->this].step != fields.field[i].step + 1) {
                goto BREAK;
            }
            parent = parent->next;
        }
        while (child) {
            if (fields.field[child->this].step != fields.field[i].step) {
                goto BREAK;
            }
            child = child->next;
        }
        while (p) {
            if (p->this == i)
                goto BREAK;
            p = p->next;
        }
        *todo = push (i, *todo);
    BREAK:;
    }
    return *todo;
}

void print_field_list (const char *title, Node *head, Fields fields) {
    PetscPrintf (PETSC_COMM_WORLD, "%s: ", title);
    while (head) {
        PetscPrintf (
            PETSC_COMM_WORLD, "%s[%zu]%s", fields.field[head->this].name,
            fields.field[head->this].step + 1, (head->next ? ", " : ""));
        head = head->next;
    }
    PetscPrintf (PETSC_COMM_WORLD, "\n\n");
}

void update (FIELD id, Fields *fields) {
    Field *field = &fields->field[id];
    field->step++;
    if (field->ncid_in)
        read_field (id, fields);
    if (field->update)
        field->update (id, fields);
}
