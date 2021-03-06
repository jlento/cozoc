#pragma once

#include <stdlib.h>

// figlet OZO | sed 's/^/"/;s/\\/\\\\/g;s/$/\\n"/'
#define BANNER                                                                 \
    "  ___ ________  \n"                                                       \
    " / _ \\__  / _ \\ \n"                                                     \
    "| | | |/ / | | |\n"                                                       \
    "| |_| / /| |_| |\n"                                                       \
    " \\___/____\\___/ \n"                                                     \
    "                \n"
#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

#if GCC_VERSION > 40899
#define let __auto_type const
#define var __auto_type
#endif


#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

#define ERROR(msg) SETERRQ (PETSC_COMM_WORLD, 1, msg)
#define info(...) PetscPrintf (PETSC_COMM_WORLD, __VA_ARGS__)

static inline size_t max_of_size_t (size_t i, size_t j) {
    return (i > j ? i : j);
}
static inline size_t min_of_size_t (size_t i, size_t j) {
    return (i < j ? i : j);
}
