#pragma once

#include "fields.h"
#include "operators.h"

#define NUM_TARGET_TYPE 2
#define NUM_TARGET 10

#define new_target_list(...)                                                \
    _new_target_list (                                                             \
        sizeof ((TARGET[]){__VA_ARGS__}) / sizeof (TARGET),                    \
        (TARGET[]){__VA_ARGS__})

enum TARGET_TYPE { TARGET_TYPE_FIELD, TARGET_TYPE_OPERATOR };

enum TARGET {
    TARGET_FIELD_DIABATIC_HEATING,        // Context
    TARGET_FIELD_H_DIABATIC,              // Input
    TARGET_FIELD_OMEGA_Q,                 // Result
    TARGET_FIELD_RTHBLTEN,                // Input
    TARGET_FIELD_RTHCUTEN,                // Input
    TARGET_FIELD_RTHRATEN,                // Input
    TARGET_FIELD_TEMPERATURE,             // Input
    TARGET_FIELD_WIND_U,                  // Input
    TARGET_FIELD_WIND_V,                  // Input
    TARGET_OPERATOR_GENERALIZED_OMEGA,    // Context
};

typedef enum TARGET_TYPE TARGET_TYPE;
typedef enum TARGET      TARGET;
typedef struct Targets   Targets;
typedef struct Target    Target;
typedef struct TARGETS   TARGETS;

struct Target {
    TARGET_TYPE type;
    union {
        Field field;
        Operator operator;
    } target;
    size_t time;
};

struct Targets {
    Context context;
    Target target[NUM_TARGET];
};

struct TARGETS {
    TARGET this;
    TARGETS *next;
};

Targets new_targets (Options options, NCFile ncfile);
TARGETS *push (TARGET target, TARGETS *oldhead);
TARGETS *pop (TARGETS **head);
TARGETS *_new_target_list (const size_t, const TARGET[]);
