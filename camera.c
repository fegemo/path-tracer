#include "camera.h"
#include "geom.h"
#include <math.h>

///
/// updates the camera orthogonal basis according to its origin (orig) and target
///
void updateCameraBasis(Camera *camera) {
	float xDirection = sin(camera->yaw) * cos(camera->pitch);
	float yDirection = sin(camera->pitch);
	float zDirection = cos(camera->yaw) * cos(camera->pitch);

	const vec directionToCamera = {xDirection, yDirection, zDirection};
	vsmul(camera->dir, -1, directionToCamera);
	vnorm(camera->dir);

	const vec up = {0.f, 1.f, 0.f};
	const float fov = (FLOAT_PI / 180.f) * HARD_CODED_CAMERA_FOV;

	// axis x is orthogonal to the camera direction and up
	vxcross(camera->x, camera->dir, up);
	vnorm(camera->x);
	// multiplies x axis by aspectRatio * fov
	vsmul(camera->x, camera->width * fov / camera->height, camera->x);

	// axis y is orthogonal to x and camera direction
	vxcross(camera->y, camera->x, camera->dir);
	vnorm(camera->y);

	// multiplies y axis by the fov
	vsmul(camera->y, fov, camera->y);
}


