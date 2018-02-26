#pragma once

#include "context.h"
#include "defs.h"
#include "io.h"
#include "options.h"
#include <netcdf.h>
#include <petscvec.h>

typedef enum FIELD {
    FIELD_RTHCUTEN,            // Input
    FIELD_RTHRATEN,            // Input
    FIELD_RTHBLTEN,            // Input
    FIELD_H_DIABATIC,          // Input
    FIELD_DIABATIC_HEATING,    // Context
    FIELD_OMEGA_Q              // Result
} FIELD;
#define NUM_FIELD 6

typedef struct Field  Field;
typedef struct Fields Fields;
typedef struct Node   Node;
typedef void (*UpdateFieldFn) (FIELD, Fields *);

struct Field {
    size_t        step;
    int           ncid_in;
    int           ncid_out;
    char          name[NC_MAX_NAME + 1];
    char *        description;
    char *        units;
    Vec           vec;
    Node *        parents;
    Node *        children;
    UpdateFieldFn update;
};

struct Fields {
    Context ctx;
    Field   field[NUM_FIELD];
};

struct Node {
    FIELD this;
    Node *next;
};

Fields new_fields (Options, NCFile);
Node *more_todo (Fields, Node **);

Node *push (FIELD, Node *);
Node *pop (Node **);
void  for_all (Node *, void (*callback) (FIELD));
void  draw_tree (Fields, const char *);
void print_field_list (const char *title, Node *head, Fields fields);

#define new_list(...)                                                          \
    ({                                                                         \
        Node * p       = 0;                                                    \
        FIELD  arr[]   = {__VA_ARGS__};                                        \
        size_t num_arr = ARRAY_SIZE (arr);                                     \
        for (size_t i = 0; i < num_arr; i++) {                                 \
            p = push (arr[i], p);                                              \
        }                                                                      \
        p;                                                                     \
    })

void update (FIELD, Fields *);
