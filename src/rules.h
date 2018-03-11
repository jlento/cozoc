#pragma once

#include "context.h"
#include "targets.h"

typedef struct Rule  Rule;
typedef struct Rules Rules;

struct Rule {
    TARGET  target;
    TARGETS *prerequisites;
    void (*recipe) (Context *);
};

struct Rules {
    Rule rule[NUM_TARGET];
};

Rules new_rules (void);
void  run (const Rules*, Targets*);
void  draw (const Rules*, const Targets*, const char *);
