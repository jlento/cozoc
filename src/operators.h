#pragma once

#include "netcdf.h"
#include <petscksp.h>

struct Operator {
    const char  name[NC_MAX_NAME + 1];
    const char *description;
    Mat         mat;
};

typedef struct Operator Operator;
