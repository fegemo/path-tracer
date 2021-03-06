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
#include "scene.h"

int workGroupSize = 1;

static vec *colors;
static unsigned int *seeds;
Camera camera;
int currentSample = 0;
Object *objects;
Material *materials;
unsigned int objectCount;
unsigned int lightCount;
Scene scene;
double startRenderingTime;

extern char captionLine1[];
extern char captionLine2[];

void updateRenderingStatistics(double);


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
	double frameStartTime = wallClockTime();

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
			vec radiance;
			RadiancePathTracing(objects, objectCount, lightCount, &ray,
					&seeds[i2], &seeds[i2 + 1], scene.skyColor1, scene.skyColor2, materials, &radiance);
//			RadianceDirectLighting(objects, objectCount, lightCount, &ray,
//					&seeds[i2], &seeds[i2 + 1], scene.skyColor1, scene.skyColor2, materials, &radiance);

			if (currentSample == 0) {
				colors[i] = radiance;
			}
			else {
				const float k1 = currentSample;
				const float k2 = 1.f / (k1 + 1.f);
				colors[i].x = (colors[i].x * k1 + radiance.x) * k2;
				colors[i].y = (colors[i].y * k1 + radiance.y) * k2;
				colors[i].z = (colors[i].z * k1 + radiance.z) * k2;
			}

			pixels[y * width + x] =
                    (toInt(colors[i].x, scene.gammaCorrection)) |
					(toInt(colors[i].y, scene.gammaCorrection) << 8) |
					(toInt(colors[i].z, scene.gammaCorrection) << 16);
		}
	}

	updateRenderingStatistics(frameStartTime);

	currentSample++;
}

/// Calculates the current rendering statistics: ellapsed time (beginning and this frame),
/// total passes, samples per second
///
/// Rendering time: total (this frame)
/// Passes: total (per second)
void updateRenderingStatistics(double frameStartTime) {
	// time spent on this "glut frame"... updates the string samples/sec to render it
	const double elapsedTimeThisFrame = wallClockTime() - frameStartTime;
	const double elapsedTime = wallClockTime() - startRenderingTime;
	char timeSinceBeginning[30];
	getHumanReadableTime(elapsedTime, timeSinceBeginning);
	sprintf(captionLine1, "Time:   %s  (%.3fs frames/s)", timeSinceBeginning, elapsedTimeThisFrame);

	const double sampleSec = height * width / elapsedTimeThisFrame;
	sprintf(captionLine2, "Passes: %d  (%.1fK samples/s)", currentSample, sampleSec / 1000.f);
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

	updateCameraBasis(&camera);
	currentSample = 0;
	startRenderingTime = wallClockTime();
	updateRendering();
}

int main(int argc, char *argv[]) {
	fprintf(stderr, "Usage: %s\n", argv[0]);
	fprintf(stderr, "Usage: %s <window width> <window height> <scene file>\n", argv[0]);

    char *sceneName;
	if (argc == 4) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
		sceneName = argv[3];
	} else if (argc == 1) {
		width = 480;
		height = 320;
		sceneName = "scenes/cornell-original.txt";
	} else {
		exit(-1);
	}
    readScene(sceneName);

	updateCameraBasis(&camera);

	printf("About to allocate opencl buffers...\n");
	allocateBuffers();

	printf("About to init glut...\n");
	char windowTitle[150];
	sprintf(windowTitle, "FegemoPT: %s (%d objects)", sceneName, objectCount);
	initGlut(argc, argv, windowTitle);

	printf("About to start the main loop...\n");
    glutMainLoop( );

	return 0;
}
