#ifndef WRFNC_H
#define WRFNC_H

enum dimensions {TIME,ZDIM,YDIM,XDIM,NDIMS};
extern const char *dimnames[NDIMS];

enum fields {z,F,NFIELDS};
extern const char *fieldnames[NFIELDS];

#endif /* WRFNC_H */
