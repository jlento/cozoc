#pragma once

#include "context.h"
#include "targets.h"

typedef struct Rule  Rule;
typedef struct Rules Rules;

struct Rule {
    TARGETS *prerequisites;
    void (*recipe) (TARGET, Targets *, const Rules *, Context *);
};

struct Rules {
    Rule rule[NUM_TARGET];
};

Rules new_rules (void);
void  run (const Rules *, Targets *, Context *);
void  draw (const Rules *, const Targets *, const char *);
