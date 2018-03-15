#pragma once

#include "context.h"
#include "defs.h"
#include "io.h"
#include "options.h"
#include <netcdf.h>
#include <petscvec.h>

enum FIELD_TYPE {
    FIELD_TYPE_INPUT,
    FIELD_TYPE_INTERMEDIATE,
    FIELD_TYPE_OUTPUT
};

typedef enum FIELD_TYPE FIELD_TYPE;

struct Field {
    int         ncid;
    const char  name[NC_MAX_NAME + 1];
    const char *description;
    const char *units;
    Vec         vec;
};

typedef struct Field Field;
