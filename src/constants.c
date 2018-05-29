#include "constants.h"
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

const double Gravitational_acceleration       = 9.80665;    // (m/s^2)
const double Specific_gas_constant_of_dry_air = 287.058;    // R
const double Specific_heat_of_dry_air         = 1004.0;     // c_p
const double sigmamin                         = 2.0e-7;
const double etamin                           = 2.0e-6;

#define SIZE_BTBUF 100
void print_backtrace (void) {
    int j, nptrs;
    void * buffer[SIZE_BTBUF];
    char **strings;

    nptrs = backtrace (buffer, SIZE_BTBUF);
    printf ("backtrace() returned %d addresses\n", nptrs);

    /* The call backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO)
       would produce similar output to the following: */

    strings = backtrace_symbols (buffer, nptrs);

    if (strings == NULL) {
        perror ("backtrace_symbols");
        exit (EXIT_FAILURE); }

    for (j = 0; j < nptrs; j++)
        printf ("%s\n", strings[j]);

    free (strings); }
