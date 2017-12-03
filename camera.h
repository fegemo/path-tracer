#ifndef _CAMERA_H
#define	_CAMERA_H

#include "vec.h"
#define HARD_CODED_CAMERA_FOV 45.f

typedef struct {
	// values defined by the scene/user interaction
	vec orig, target;
	// values used to facilitate interaction
	float pitch, yaw;
	// resolution in pixels
	float width, height;

	// values calculated from orig/target (direction, x and y basis vectors)
	vec dir, x, y;

} Camera;

void updateCameraBasis(Camera *);

#endif	/* _CAMERA_H */

