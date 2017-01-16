#ifndef WRFNC_H
#define WRFNC_H

enum dimensions { TIME, ZDIM, YDIM, XDIM, NDIMS };
extern const char* dimnames[NDIMS];

enum fields { TIME_COORDINATE, Z_COORDINATE, FRICTION, NFIELDS };
extern const char* fieldnames[NFIELDS];

#endif /* WRFNC_H */
