/*
Copyright (c) 2009 David Bucciarelli (davibu@interfree.it)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
 * Based on smallpt, a Path Tracer by Kevin Beason, 2008
 * Modified by David Bucciarelli to show the output via OpenGL/GLUT, ported
 * to C, work with float, fixed RR, ported to OpenCL, etc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "camera.h"
#include "simplernd.h"
#include "geom.h"
#include "geomfunc.h"
#include "displayfunc.h"

int workGroupSize = 1;

static vec *colors;
static unsigned int *seeds;
Camera camera;
static int currentSample = 0;
Object *objects;
unsigned int objectCount;

void freeBuffers() {
	free(seeds);
	free(colors);
	free(pixels);
}

void allocateBuffers() {
	const int pixelCount = height * width ;
	int i;
	colors = malloc(sizeof(vec[pixelCount]));

	seeds = malloc(sizeof(unsigned int[pixelCount * 2]));
	for (i = 0; i < pixelCount * 2; i++) {
		seeds[i] = rand();
		if (seeds[i] < 2)
			seeds[i] = 2;
	}

	pixels = malloc(sizeof(unsigned int[pixelCount]));
}

void updateRendering(void) {
	double startTime = wallClockTime();

	const float invWidth = 1.f / width;
	const float invHeight = 1.f / height;

	int x, y;
	for (y = 0; y < height; y++) { /* Loop over image rows */
		for (x = 0; x < width; x++) { /* Loop cols */
			const int i = (height - y - 1) * width + x;
			const int i2 = 2 * i;

			const float r1 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
			const float r2 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
			const float kcx = (x + r1) * invWidth - .5f;
			const float kcy = (y + r2) * invHeight - .5f;

			vec rdir;
			vinit(rdir,
				camera.x.x * kcx + camera.y.x * kcy + camera.dir.x,
				camera.x.y * kcx + camera.y.y * kcy + camera.dir.y,
				camera.x.z * kcx + camera.y.z * kcy + camera.dir.z);

			vec rorig;
			vsmul(rorig, 0.1f, rdir);
			vadd(rorig, rorig, camera.orig)

			vnorm(rdir);
			const Ray ray = {rorig, rdir};
			vec r;
			RadiancePathTracing(objects, objectCount, &ray,
					&seeds[i2], &seeds[i2 + 1], &r);

			if (currentSample == 0)
				colors[i] = r;
			else {
				const float k1 = currentSample;
				const float k2 = 1.f / (k1 + 1.f);
				colors[i].x = (colors[i].x * k1 + r.x) * k2;
				colors[i].y = (colors[i].y * k1 + r.y) * k2;
				colors[i].z = (colors[i].z * k1 + r.z) * k2;
			}

			pixels[y * width + x] = toInt(colors[i].x) |
					(toInt(colors[i].y) << 8) |
					(toInt(colors[i].z) << 16);
		}
	}

	const float elapsedTime = wallClockTime() - startTime;
	const float sampleSec = height * width / elapsedTime;
	sprintf(captionBuffer, "Rendering time %.3f sec (pass %d)  Sample/sec  %.1fK\n",
		elapsedTime, currentSample, sampleSec / 1000.f);

	currentSample++;
}

void reInitSceneObjects() {
	currentSample = 0;
}

void reInitViewpointDependentBuffers(const int reallocBuffers) {
	// Check if I have to reallocate buffers
	if (reallocBuffers) {
		freeBuffers();
		allocateBuffers();
	}

	updateCamera();
	currentSample = 0;
	updateRendering();
}

int main(int argc, char *argv[]) {
	amiSmallptCPU = 1;

	fprintf(stderr, "Usage: %s\n", argv[0]);
	fprintf(stderr, "Usage: %s <window width> <window height> <scene file>\n", argv[0]);

	if (argc == 4) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
		readScene(argv[3]);
	} else if (argc == 1) {
		width = 480;
		height = 320;
		readScene("scenes/cornell.scn");
	} else {
		exit(-1);
	}

	updateCamera();

	/*------------------------------------------------------------------------*/

	allocateBuffers();

	/*------------------------------------------------------------------------*/

	initGlut(argc, argv, "SmallPT CPU V1.6 (Written by David Bucciarelli)");

    glutMainLoop( );

	return 0;
}
