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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#if defined(__linux__) || defined(__APPLE__)
#include <sys/time.h>
#elif defined (WIN32)
#include <windows.h>
#else
        Unsupported Platform !!!
#endif

#include "camera.h"
#include "geom.h"
#include "displayfunc.h"
#define TINYOBJ_LOADER_C_IMPLEMENTATION
#include "tinyobj_loader_c.h"

extern void reInitViewpointDependentBuffers(const int);
extern void reInitSceneObjects();
extern void updateRendering();
extern void updateCamera();

extern Camera camera;
extern Object *objects;
extern unsigned int objectCount;

int amiSmallptCPU;

int width = 640;
int height = 480;
unsigned int *pixels;
char captionBuffer[256];

static int shouldRenderHelp = 1;
//static int currentSphere;

double wallClockTime() {
#if defined(__linux__) || defined(__APPLE__)
	struct timeval t;
	gettimeofday(&t, NULL);

	return t.tv_sec + t.tv_usec / 1000000.0;
#elif defined (WIN32)
	return GetTickCount() / 1000.0;
#else
	Unsupported Platform !!!
#endif
}

static const char* readFile(size_t* fileLength, const char* fileName) {
    char *contents;
    FILE *file = fopen(fileName, "rb");
    fseek(file, 0, SEEK_END);
    *fileLength = ftell(file);
    rewind(file);
    contents = malloc(*fileLength * (sizeof(char)));
    fread(contents, sizeof(char), *fileLength, file);
    fclose(file);

    return contents;
}

///
/// Renders text to the opengl window
///
static void renderString(void *font, const char *string) {
	int len, i;

	len = (int)strlen(string);
	for (i = 0; i < len; i++) {
		glutBitmapCharacter(font, string[i]);
	}
}

///
/// Renders the help box to the opengl window
///
static void renderHelp() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4f(0.f, 0.f, 0.5f, 0.5f);
	glRecti(40, 40, 600, 440);

	glColor3f(1.f, 1.f, 1.f);
	glRasterPos2i(300, 420);
	renderString(GLUT_BITMAP_HELVETICA_18, "Help");

	glRasterPos2i(60, 390);
	renderString(GLUT_BITMAP_HELVETICA_18, "h - toggle Help");
	glRasterPos2i(60, 360);
	renderString(GLUT_BITMAP_HELVETICA_18, "arrow Keys - rotate camera left/right/up/down");
	glRasterPos2i(60, 330);
	renderString(GLUT_BITMAP_HELVETICA_18, "a and d - move camera left and right");
	glRasterPos2i(60, 300);
	renderString(GLUT_BITMAP_HELVETICA_18, "w and s - move camera forward and backward");
	glRasterPos2i(60, 270);
	renderString(GLUT_BITMAP_HELVETICA_18, "r and f - move camera up and down");
	glRasterPos2i(60, 240);
	renderString(GLUT_BITMAP_HELVETICA_18, "PageUp and PageDown - move camera target up and down");
	glRasterPos2i(60, 210);
	renderString(GLUT_BITMAP_HELVETICA_18, "+ and - - to select next/previous object");
	glRasterPos2i(60, 180);
	renderString(GLUT_BITMAP_HELVETICA_18, "2, 3, 4, 5, 6, 8, 9 - to move selected object");

	glDisable(GL_BLEND);
}

///
/// Reads a scene file and saves the objects in the objects global array.
///
void readScene(char *fileName) {
	fprintf(stderr, "Reading scene: %s\n", fileName);

	FILE *f = fopen(fileName, "r");
	if (!f) {
		fprintf(stderr, "Failed to open file: %s\n", fileName);
		exit(-1);
	}

	// line 1: camera configuration: camera eye.x eye.y eye.z  target.x target.y target.z
	int c = fscanf(f,"camera %f %f %f  %f %f %f\n",
			&camera.orig.x, &camera.orig.y, &camera.orig.z,
			&camera.target.x, &camera.target.y, &camera.target.z);
	if (c != 6) {
		fprintf(stderr, "Failed to read 6 camera parameters: %d\n", c);
		exit(-1);
	}

	// line 2: object counts: size n
	c = fscanf(f,"size %u\n", &objectCount);
	if (c < 1) {
		fprintf(stderr, "Failed to read object count: %d\n", c);
		exit(-1);
	}
	fprintf(stderr, "Scene size: %d\n", objectCount);

	// line 3+: object descriptors
	printf("Allocated objects to have %d positions.\n", objectCount);
	objects = (Object *)malloc(sizeof(Object) * objectCount);
	unsigned int lineIndex;
	unsigned int addedObjects = 0;
	for (lineIndex = 0; lineIndex < objectCount; lineIndex++) {
        Object *obj = &objects[lineIndex + addedObjects];
		int mat;
		char objectTypeString[100];

		int readValuesCount = fscanf(f, "%s", objectTypeString);
        if (readValuesCount != 1) {
            fprintf(stderr, "Failed to read object description type. Expected 1, found %d value(s)\n", readValuesCount);
            exit(-1);
        }
        printf("reading an object of type %s and saving on position %d.\n", objectTypeString, lineIndex + addedObjects);

		if (strcmp(objectTypeString, "sphere") == 0) {
		    obj->type = SPHERE;
            // line x: sphere     radius  xyz  emission  color  reflection
            readValuesCount = fscanf(f,"%f  %f %f %f  %f %f %f  %f %f %f  %d\n",
				&obj->radius,
				&obj->center.x, &obj->center.y, &obj->center.z,
				&obj->emission.x, &obj->emission.y, &obj->emission.z,
				&obj->color.x, &obj->color.y, &obj->color.z,
				&mat);
            switch (mat) {
                case 0: obj->refl = DIFF; break;
                case 1: obj->refl = SPEC; break;
                case 2: obj->refl = REFR; break;
                default:
                    fprintf(stderr, "Failed to read material type for sphere #%d: %d\n", lineIndex, mat);
                    exit(-1);
                    break;
            }

            if (readValuesCount != 11) {
                fprintf(stderr, "Failed to read sphere #%d: %d\n", lineIndex, readValuesCount);
                exit(-1);
            }

		} else if (strcmp(objectTypeString, "triangle") == 0) {
		    obj->type = TRIANGLE;
            readValuesCount = fscanf(f, "%f %f %f  %f %f %f  %f %f %f  %f %f %f  %f %f %f  %d\n",
                &obj->p1.x, &obj->p1.y, &obj->p1.z,
                &obj->p2.x, &obj->p2.y, &obj->p2.z,
                &obj->p3.x, &obj->p3.y, &obj->p3.z,
                &obj->emission.x, &obj->emission.y, &obj->emission.z,
                &obj->color.x, &obj->color.y, &obj->color.z,
                &mat);
//            printf("triangle: %f %f %f   %f %f %f   %f %f %f   %f %f %f   %f %f %f    %d\n",
//                   obj->p1.x, obj->p1.y, obj->p1.z, obj->p2.x, obj->p2.y, obj->p2.z, obj->p3.x, obj->p3.y, obj->p3.z, obj->emission.x, obj->emission.y, obj->emission.z, obj->color.x, obj->color.y, obj->color.z, obj->refl);
            switch (mat) {
                case 0: obj->refl = DIFF; break;
                case 1: obj->refl = SPEC; break;
                case 2: obj->refl = REFR; break;
                    fprintf(stderr, "Failed to read material type for triangle #%d: %d\n", lineIndex, mat);
                    exit(-1);
                    break;
            }

            if (readValuesCount != 16) {
                fprintf(stderr, "Failed to read triangle #%d: %d\n", lineIndex, readValuesCount);
                exit(-1);
            }
		} else if (strcmp(objectTypeString, "model") == 0) {
		    char modelFilePath[200];
		    vec modelPosition;
            vec modelScale;
		    obj->type = MODEL;
		    readValuesCount = fscanf(f, "%s  %f %f %f  %f %f %f  %f %f %f  %f %f %f  %d\n", modelFilePath,
                               &modelPosition.x, &modelPosition.y, &modelPosition.z,
                               &modelScale.x, &modelScale.y, &modelScale.z,
                               &obj->emission.x, &obj->emission.y, &obj->emission.z,
                               &obj->color.x, &obj->color.y, &obj->color.z,
                               &mat);

            if (readValuesCount != 14) {
                fprintf(stderr, "Failed to read model #%d: %d\n", lineIndex, readValuesCount);
                exit(-1);
            }

            switch (mat) {
                case 0: obj->refl = DIFF; break;
                case 1: obj->refl = SPEC; break;
                case 2: obj->refl = REFR; break;
                    fprintf(stderr, "Failed to read material type for model #%d: %d\n", lineIndex, mat);
                    exit(-1);
                    break;
            }

            tinyobj_attrib_t attrib;
            tinyobj_shape_t* shapes = NULL;
            size_t num_shapes;
            tinyobj_material_t* materials = NULL;
            size_t num_materials;
            size_t data_len = 0;

            const char* data = readFile(&data_len, modelFilePath);
            if (data == NULL) {
                exit(-1);
            }

//            unsigned int flags = TINYOBJ_FLAG_TRIANGULATE;
            unsigned int flags = NULL;
            int ret = tinyobj_parse_obj(&attrib, &shapes, &num_shapes, &materials,
                                    &num_materials, data, data_len, flags);
            if (ret != TINYOBJ_SUCCESS) {
                fprintf(stderr, "Error loading obj. Exiting...");
                exit(1);
            }

            printf("# of shapes    = %d\n", (int)num_shapes);
            printf("# of materials = %d\n", (int)num_materials);


            unsigned int numTriangles = attrib.num_face_num_verts;
            printf("Started reading model with %d triangles.\n", numTriangles);

            objects = (Object *) realloc(objects, sizeof(Object) * (objectCount + addedObjects + numTriangles));
            printf("Reallocated objects to have %d positions.\n", objectCount + addedObjects + numTriangles);
            obj = &objects[lineIndex + addedObjects];
            unsigned int faceIndex;
            unsigned int faceOffset = 0;
            for (faceIndex = 0; faceIndex < numTriangles; faceIndex++) { //6, 7
                Object *tri = &objects[lineIndex + addedObjects + faceIndex + 1];
//                printf("adding new triangle and saving on position %d.\n", lineIndex + addedObjects + faceIndex  + 1);
                tri->type = TRIANGLE;
                tri->emission = obj->emission;
                tri->color = obj->color;
                tri->refl = obj->refl;

                int idxVertex1 = attrib.faces[faceOffset+0].v_idx;
                int idxVertex2 = attrib.faces[faceOffset+1].v_idx;
                int idxVertex3 = attrib.faces[faceOffset+2].v_idx;

                vec p1; vinit(p1, attrib.vertices[3*idxVertex1+0], attrib.vertices[3*idxVertex1+1], attrib.vertices[3*idxVertex1+2]);
                vec p2; vinit(p2, attrib.vertices[3*idxVertex2+0], attrib.vertices[3*idxVertex2+1], attrib.vertices[3*idxVertex2+2]);
                vec p3; vinit(p3, attrib.vertices[3*idxVertex3+0], attrib.vertices[3*idxVertex3+1], attrib.vertices[3*idxVertex3+2]);

                // scales the model (according to .scn file)
                vmul(p1, modelScale, p1);
                vmul(p2, modelScale, p2);
                vmul(p3, modelScale, p3);

                // translates the vertices to the model position (according to .scn file)
                vadd(p1, p1, modelPosition);
                vadd(p2, p2, modelPosition);
                vadd(p3, p3, modelPosition);

                // saves the vertices coords into the triangle structure
                tri->p1 = p1;
                tri->p2 = p2;
                tri->p3 = p3;

                // sets a center of the triangle (in case it's a light source) and its radius
                vclr(tri->center)
                vadd(tri->center, p1, tri->center);
                vadd(tri->center, p2, tri->center);
                vadd(tri->center, p3, tri->center);
                vsmul(tri->center, 1/3.0f, tri->center);

                // PAREI AQUI..... PREENCHENDO RADIUS PRA VER SE FUNCIONA LUZ RETANGULAR
                tri->radius = 7;//max(max(fabs(p1.x )))


                faceOffset += 3;
            }

            addedObjects += numTriangles;
		}

	}

	objectCount += addedObjects;

	printf("Finished parsing scene. Resulting objectCount = %d\n", objectCount);

//	int i;
	for (int i = 0; i < objectCount; i++) {
        printf("Object #%d: %d\n", i, objects[i].type);
	}

	fclose(f);
}

///
/// updates the camera orthogonal basis according to its origin (orig) and target
///
void updateCamera() {
	vsub(camera.dir, camera.target, camera.orig);
	vnorm(camera.dir);

	const vec up = {0.f, 1.f, 0.f};
	const float fov = (M_PI / 180.f) * 45.f;

	// axis x is orthogonal to the camera direction and up
	vxcross(camera.x, camera.dir, up);
	vnorm(camera.x);
	// multiplies x axis by aspectRatio * fov
	vsmul(camera.x, width * fov / height, camera.x);

	// axis y is orthogonal to x and camera direction
	vxcross(camera.y, camera.x, camera.dir);
	vnorm(camera.y);
	// multiplies y axis by the fov
	vsmul(camera.y, fov, camera.y);
}

///
/// updates the opengl window by asking opencl to execute the kernel and then
/// asking glut to redraw the window
///
void update(void) {
	updateRendering();
	glutPostRedisplay();
}

///
/// draws the window by looking at the color buffer that was "outputted" from opencl
///
void displayScene(void) {
    // erases the screen and draws the color buffer (pixels array)
	glClear(GL_COLOR_BUFFER_BIT);
	glRasterPos2i(0, 0);
	glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

	// renders the title at the top
	glColor3f(1.f, 1.f, 1.f);
	glRasterPos2i(4, height - 16);
//	if (amiSmallptCPU)
//		renderString(GLUT_BITMAP_HELVETICA_18, "SmallptCPU v1.6 (Written by David Bucciarelli)");
//	else
//		renderString(GLUT_BITMAP_HELVETICA_18, "SmallptGPU v1.6 (Written by David Bucciarelli)");

	// renders the stats at the bottom
	glColor3f(1.f, 1.f, 1.f);
	glRasterPos2i(4, 10);
	renderString(GLUT_BITMAP_HELVETICA_18, captionBuffer);

	// renders the help box with its instructions
//	if (shouldRenderHelp) {
//		renderHelp();
//	}

	glutSwapBuffers();
}

///
/// Resizes the width/height of the image and recalculates the camera.
///
void reshape(int newWidth, int newHeight) {
	width = newWidth;
	height = newHeight;

	glViewport(0, 0, width, height);
	glLoadIdentity();
	glOrtho(0.f, width - 1.f, 0.f, height - 1.f, -1.f, 1.f);

	reInitViewpointDependentBuffers(1);

	glutPostRedisplay();
}

#define MOVE_STEP 10.0f
#define ROTATE_STEP (2.f * M_PI / 180.f)
void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 'p': {
			FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
			if (!f) {
				fprintf(stderr, "Failed to open image file: image.ppm\n");
			} else {
				fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);

				int x, y;
				for (y = height - 1; y >= 0; --y) {
					unsigned char *p = (unsigned char *)(&pixels[y * width]);
					for (x = 0; x < width; ++x, p += 4)
						fprintf(f, "%d %d %d ", p[0], p[1], p[2]);
				}

				fclose(f);
			}
			break;
		}
		case 27: /* Escape key */
			fprintf(stderr, "Done.\n");
			exit(0);
			break;
		case ' ': /* Refresh display */
			reInitViewpointDependentBuffers(1);
			break;
		case 'a': {
			vec dir = camera.x;
			vnorm(dir);
			vsmul(dir, -MOVE_STEP, dir);
			vadd(camera.orig, camera.orig, dir);
			vadd(camera.target, camera.target, dir);
			reInitViewpointDependentBuffers(0);
			break;
		}
		case 'd': {
			vec dir = camera.x;
			vnorm(dir);
			vsmul(dir, MOVE_STEP, dir);
			vadd(camera.orig, camera.orig, dir);
			vadd(camera.target, camera.target, dir);
			reInitViewpointDependentBuffers(0);
			break;
		}
		case 'w': {
			vec dir = camera.dir;
			vsmul(dir, MOVE_STEP, dir);
			vadd(camera.orig, camera.orig, dir);
			vadd(camera.target, camera.target, dir);
			reInitViewpointDependentBuffers(0);
			break;
		}
		case 's': {
			vec dir = camera.dir;
			vsmul(dir, -MOVE_STEP, dir);
			vadd(camera.orig, camera.orig, dir);
			vadd(camera.target, camera.target, dir);
			reInitViewpointDependentBuffers(0);
			break;
		}
		case 'r':
			camera.orig.y += MOVE_STEP;
			camera.target.y += MOVE_STEP;
			reInitViewpointDependentBuffers(0);
			break;
		case 'f':
			camera.orig.y -= MOVE_STEP;
			camera.target.y -= MOVE_STEP;
			reInitViewpointDependentBuffers(0);
			break;
//		case '+':
//			currentSphere = (currentSphere + 1) % sphereCount;
//			fprintf(stderr, "Selected sphere %d (%f %f %f)\n", currentSphere,
//					spheres[currentSphere].p.x, spheres[currentSphere].p.y, spheres[currentSphere].p.z);
//			reInitSceneObjects();
//			break;
//		case '-':
//			currentSphere = (currentSphere + (sphereCount - 1)) % sphereCount;
//			fprintf(stderr, "Selected sphere %d (%f %f %f)\n", currentSphere,
//					spheres[currentSphere].p.x, spheres[currentSphere].p.y, spheres[currentSphere].p.z);
//			reInitSceneObjects();
//			break;
//		case '4':
//			spheres[currentSphere].p.x -= 0.5f * MOVE_STEP;
//			reInitSceneObjects();
//			break;
//		case '6':
//			spheres[currentSphere].p.x += 0.5f * MOVE_STEP;
//			reInitSceneObjects();
//			break;
//		case '8':
//			spheres[currentSphere].p.z -= 0.5f * MOVE_STEP;
//			reInitSceneObjects();
//			break;
//		case '2':
//			spheres[currentSphere].p.z += 0.5f * MOVE_STEP;
//			reInitSceneObjects();
//			break;
//		case '9':
//			spheres[currentSphere].p.y += 0.5f * MOVE_STEP;
//			reInitSceneObjects();
//			break;
//		case '3':
//			spheres[currentSphere].p.y -= 0.5f * MOVE_STEP;
//			reInitSceneObjects();
//			break;
		case 'h':
			shouldRenderHelp = (!shouldRenderHelp);
			break;
		default:
			break;
	}
}

void specialKeyboard(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_UP: {
			vec t = camera.target;
			vsub(t, t, camera.orig);
			t.y = t.y * cos(-ROTATE_STEP) + t.z * sin(-ROTATE_STEP);
			t.z = -t.y * sin(-ROTATE_STEP) + t.z * cos(-ROTATE_STEP);
			vadd(t, t, camera.orig);
			camera.target = t;
			reInitViewpointDependentBuffers(0);
			break;
		}
		case GLUT_KEY_DOWN: {
			vec t = camera.target;
			vsub(t, t, camera.orig);
			t.y = t.y * cos(ROTATE_STEP) + t.z * sin(ROTATE_STEP);
			t.z = -t.y * sin(ROTATE_STEP) + t.z * cos(ROTATE_STEP);
			vadd(t, t, camera.orig);
			camera.target = t;
			reInitViewpointDependentBuffers(0);
			break;
		}
		case GLUT_KEY_LEFT: {
			vec t = camera.target;
			vsub(t, t, camera.orig);
			t.x = t.x * cos(-ROTATE_STEP) - t.z * sin(-ROTATE_STEP);
			t.z = t.x * sin(-ROTATE_STEP) + t.z * cos(-ROTATE_STEP);
			vadd(t, t, camera.orig);
			camera.target = t;
			reInitViewpointDependentBuffers(0);
			break;
		}
		case GLUT_KEY_RIGHT: {
			vec t = camera.target;
			vsub(t, t, camera.orig);
			t.x = t.x * cos(ROTATE_STEP) - t.z * sin(ROTATE_STEP);
			t.z = t.x * sin(ROTATE_STEP) + t.z * cos(ROTATE_STEP);
			vadd(t, t, camera.orig);
			camera.target = t;
			reInitViewpointDependentBuffers(0);
			break;
		}
		case GLUT_KEY_PAGE_UP:
			camera.target.y += MOVE_STEP;
			reInitViewpointDependentBuffers(0);
			break;
		case GLUT_KEY_PAGE_DOWN:
			camera.target.y -= MOVE_STEP;
			reInitViewpointDependentBuffers(0);
			break;
		default:
			break;
	}
}

///
/// Initializes glut and opengl to sensible defaults.
///
void initGlut(int argc, char *argv[], char *windowTittle) {
    printf("About to init window stuff...\n");
    glutInitWindowSize(width, height);
    glutInitWindowPosition(0,0);
    printf("About to init glut displaymode...\n");
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    printf("About to init glut...\n");
	glutInit(&argc, argv);

    printf("About to create window...\n");
	glutCreateWindow(windowTittle);

    printf("About to register callbacks...\n");
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(specialKeyboard);
    glutDisplayFunc(displayScene);
	glutIdleFunc(update);

    printf("About to config viewport and projection...\n");
	glViewport(0, 0, width, height);
	glLoadIdentity();
	glOrtho(0.f, width - 1.f, 0.f, height - 1.f, -1.f, 1.f);
}
