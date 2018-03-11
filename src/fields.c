#include "fields.h"
#include "context.h"
#include "defs.h"
#include <petscsys.h>
#include <stdio.h>

/*
void read_field (FIELD id, Fields *fields) {
    Field field = fields->field[id];
    file_read_3d (field.ncid_in, field.step, field.name, field.vec);
}
*/

/*
void update (FIELD id, Fields *fields) {
    Field *field = &fields->field[id];
    field->step++;
    if (field->ncid_in)
        read_field (id, fields);
    if (field->update)
        field->update (id, fields);
}
*/
