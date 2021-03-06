#ifndef _GEOM_H
#define	_GEOM_H

#include "vec.h"

#define EPSILON 0.01f
#define FLOAT_PI 3.14159265358979323846f
#define PI_OVER_TWO 1.57079632679489661923f
#define TWO_PI 6.28318530717958647693f
#define FOUR_PI 12.56637061435917295385f
#define INV_PI 0.31830988618379067154f
#define INV_TWOPI 0.15915494309189533577f

typedef struct {
	vec o, d;
} Ray;

#define rinit(r, a, b) { vassign((r).o, a); vassign((r).d, b); }
#define rassign(a, b) { vassign((a).o, (b).o); vassign((a).d, (b).d); }


enum ObjectType {
    SPHERE, TRIANGLE, MODEL
};

typedef struct {
    // for any type: type and material
    short type;
    short materialId;
    vec emission;
    vec color;

    // for spheres:
    float radius;
    vec center;

    // for triangles:
    vec p1, p2, p3;

    // for objects that are lights (ie, emit color):
    float area;

} Object;

#endif	/* _GEOM_H */

