#define GPU_KERNEL

#include "camera.h"
#include "geom.h"
#include "material.h"
#include "geomfunc.h"
#include "scene.h"

static void GenerateCameraRay(OCL_CONSTANT_BUFFER Camera *camera,
		unsigned int *seed0, unsigned int *seed1,
		// width/height of the image
		const int width, const int height,
		// (x,y) of the pixel (i,j, actually)
		const int x, const int y,
		// the output ray
		Ray *ray) {

	const float invWidth = 1.f / width;
	const float invHeight = 1.f / height;
	// gets tiny random bumps in the x,y of the pixel ([-0.5, 0.5])
	const float r1 = GetRandom(seed0, seed1) - .5f;
	const float r2 = GetRandom(seed0, seed1) - .5f;
	// finds the (x,y) of this pixel, but between [-0.5, 0.5]
	const float kcx = (x + r1) * invWidth - .5f;
	const float kcy = (y + r2) * invHeight - .5f;

	vec rdir;
	// the ray direction is the camera.dir, but tilted to point in the
	// direction of the pixel
	// for a camera with "up pointing up", a sample ray shooting at the
	// leftmost pixel, vertically centered, would be:
	vinit(rdir,
            //     -0.5*cam.x +                 0 +     cam.dir.x   = camera dir.x tilted to the left
			camera->x.x * kcx + camera->y.x * kcy + camera->dir.x,
            //              0 +                 0 +     cam.dir.y   = camera dir.y
			camera->x.y * kcx + camera->y.y * kcy + camera->dir.y,
            //              0 +                 0 +     cam.dir.z   = camera.dir.z
			camera->x.z * kcx + camera->y.z * kcy + camera->dir.z);

    // the ray origin is (direction/10) + (camera.origin)
	vec rorig;
	vsmul(rorig, 0.1f, rdir);   // 0.1 is the distance from camera to the lenses (near plane, or film)
	vadd(rorig, rorig, camera->orig)

	vnorm(rdir);
	rinit(*ray, rorig, rdir);
}

__kernel void radianceGPU(
    __global vec *colors,
    __global unsigned int *seedsInput,
	__global Object *objects,
	OCL_CONSTANT_BUFFER Camera *camera,
	const unsigned int objectCount,
	const int width, const int height,
	const int currentSample,
	__global int *pixels,
	__global int *debug,
	const unsigned int lightCount,
	const Scene scene,
	__global Material *materials) {

    debug[0] = sizeof(Object);
    const int gid = get_global_id(0);
	const int gid2 = 2 * gid;
	const int x = gid % width;
	const int y = gid / width;

	// we only need to run determine the color if the pixel lies inside the image...
	// (why would we have a y larger than the height?)
	if (y >= height)
		return;

	unsigned int seed0 = seedsInput[gid2];
	unsigned int seed1 = seedsInput[gid2 + 1];

	// generates a ray going through pixel (x,y), in a random position inside it
	// the output is returned on "ray"
	Ray ray;
	GenerateCameraRay(camera, &seed0, &seed1, width, height, x, y, &ray);

	// shoots the ray in the scene and receives back a radiance value
	// the output is returned on "radiance"
	vec radiance;
	RadiancePathTracing(objects, objectCount, lightCount, &ray, &seed0, &seed1, scene.skyColor1, scene.skyColor2, materials, &radiance);

	const int i = (height - y - 1) * width + x;
	if (currentSample == 0) {
		vassign(colors[i], radiance);
	} else {
	    // k1 holds the current iteration index (in 1 iteration, we shoot w*h rays, 1/pixel)
		const float k1 = currentSample;
		// k2 holds a number that gets smaller and smaller as the iterations increase
		const float k2 = 1.f / (currentSample + 1.f);

		// the color of the pixel is given by an weighed average of the colors of the previous
		// iterations and the current one...
		// an example:
		//  1st ray: color(1) = radiance(1)
		//  2nd ray: color(2) = (color(1) * 1 + radiance(2)) /2
		//  3rd ray: color(3) = (color(2) * 2 + radiance(3)) /3
		//  4th ray: color(4) = (color(3) * 3 + radiance(4)) /4
		//  5th ray: color(5) = (color(4) * 4 + radiance(5)) /5
		//  ...
		//  ith ray: color(i) =(color(i-1) * (i-1) + radiance(i)) /(i)
		// the average weighs the "color so far" according to how many rays we have shot for this pixel
		//
		colors[i].x = (colors[i].x * k1  + radiance.x) * k2;
		colors[i].y = (colors[i].y * k1  + radiance.y) * k2;
		colors[i].z = (colors[i].z * k1  + radiance.z) * k2;
	}

	pixels[y * width + x] =
            (toInt(colors[i].x, scene.gammaCorrection)) |         // blue
			(toInt(colors[i].y, scene.gammaCorrection) << 8) |    // green
			(toInt(colors[i].z, scene.gammaCorrection) << 16);    // red

	seedsInput[gid2] = seed0;
	seedsInput[gid2 + 1] = seed1;
}
