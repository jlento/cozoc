#pragma once

#include "context.h"
#include "defs.h"
#include "io.h"
#include "options.h"
#include <netcdf.h>
#include <petscvec.h>

struct Field {
    int         ncid_in;
    int         ncid_out;
    const char  name[NC_MAX_NAME + 1];
    const char *description;
    const char *units;
    Vec         vec;
};

typedef struct Field Field;
