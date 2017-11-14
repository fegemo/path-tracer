#ifndef _CAMERA_H
#define	_CAMERA_H

#include "vec.h"

typedef struct {
	// values defined by the scene/user interaction
	vec orig, target;
	// values calculated from orig/target (direction, x and y basis vectors)
	vec dir, x, y;
} Camera;

#endif	/* _CAMERA_H */

