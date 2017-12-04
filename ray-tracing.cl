#define GPU_KERNEL

#include "camera.h"
#include "geomfunc.h"

static void GenerateCameraRay(OCL_CONSTANT_BUFFER Camera *camera,
		unsigned int *seed0, unsigned int *seed1,
		const int width, const int height, const int x, const int y, Ray *ray) {
	const float invWidth = 1.f / width;
	const float invHeight = 1.f / height;
	const float r1 = GetRandom(seed0, seed1) - .5f;
	const float r2 = GetRandom(seed0, seed1) - .5f;
	const float kcx = (x + r1) * invWidth - .5f;
	const float kcy = (y + r2) * invHeight - .5f;

	vec rdir;
	vinit(rdir,
			camera->x.x * kcx + camera->y.x * kcy + camera->dir.x,
			camera->x.y * kcx + camera->y.y * kcy + camera->dir.y,
			camera->x.z * kcx + camera->y.z * kcy + camera->dir.z);

	vec rorig;
	vsmul(rorig, 0.1f, rdir);
	vadd(rorig, rorig, camera->orig)

	vnorm(rdir);
	rinit(*ray, rorig, rdir);
}

__kernel void radianceGPU(
    __global vec *colors, __global unsigned int *seedsInput,
	__global Object *object, OCL_CONSTANT_BUFFER Camera *camera,
	const unsigned int objectCount,
	const unsigned int lightCount,
	const int width, const int height,
	const int currentSample,
	__global int *pixels,
	__global int *debug) {

    const int gid = get_global_id(0);
	const int gid2 = 2 * gid;
	const int x = gid % width;
	const int y = gid / width;

	/* Check if we have to do something */
	if (y >= height)
		return;

	unsigned int seed0 = seedsInput[gid2];
	unsigned int seed1 = seedsInput[gid2 + 1];

	Ray ray;
	GenerateCameraRay(camera, &seed0, &seed1, width, height, x, y, &ray);

	vec radiance;
	RadianceDirectLighting(object, objectCount, lightCount, &ray, &seed0, &seed1, &radiance);

	const int i = (height - y - 1) * width + x;
	if (currentSample == 0) {
		// Jens's patch for MacOS
		vassign(colors[i], radiance);
	} else {
		const float k1 = currentSample;
		const float k2 = 1.f / (currentSample + 1.f);
		colors[i].x = (colors[i].x * k1  + radiance.x) * k2;
		colors[i].y = (colors[i].y * k1  + radiance.y) * k2;
		colors[i].z = (colors[i].z * k1  + radiance.z) * k2;
	}

	pixels[y * width + x] = toInt(colors[i].x) |
			(toInt(colors[i].y) << 8) |
			(toInt(colors[i].z) << 16);

	seedsInput[gid2] = seed0;
	seedsInput[gid2 + 1] = seed1;
}
