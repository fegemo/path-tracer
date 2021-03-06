#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <libgen.h>

#include "camera.h"
#include "geom.h"
#include "material.h"
#include "displayfunc.h"
#include "time-utils.h"
#include "scene.h"


#define TINYOBJ_LOADER_C_IMPLEMENTATION
#include "tinyobj_loader_c.h"

extern void reInitViewpointDependentBuffers(const int);
extern void reInitSceneObjects();
extern void updateRendering();

extern Camera camera;
extern Object *objects;
extern unsigned int objectCount;
extern unsigned int lightCount;
extern Material *materials;
extern Scene scene;
extern int currentSample;
extern double startRenderingTime;

char sceneTitle[100];
int width = 640;
int height = 480;
unsigned int *pixels;
char captionLine1[256];
char captionLine2[256];


// interaction with the camera
static int mouseX = 0, mouseY = 0;
static int mouseButton = 0;
static bool keyStates[256];
static bool shiftPressed = false;
static Camera originalCamera;
static int timeBefore = 0;

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

/// Removes the extension of a file name
void stripExtension(char *fileName) {
    char *end = fileName + strlen(fileName);

    while (end > fileName && *end != '.' && *end != '\\' && *end != '/') {
        --end;
    }

    if (end > fileName && *end == '.') {
        *end = '\0';
    }
}

int getMaterialType(unsigned int id, int lineIndex) {
    switch (id) {
        case 0: return LAMBERTIAN;
        case 1: return CONDUCTOR;
        case 2: return DIELECTRIC;
        case 3: return MATTE;
        case 4: return PLASTIC;
        case 5: return METAL;
        default:
            fprintf(stderr, "Failed to read material type for sphere #%d: %d\n", lineIndex, id);
            exit(-1);
    }
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

	strcpy(sceneTitle, basename(fileName));
	stripExtension(sceneTitle);

	int c;
	// part 1: gamma correction, sky colors
	c = fscanf(f, "gamma %f  sky1 %f %f %f  sky2 %f %f %f\n", &scene.gammaCorrection,
            &scene.skyColor1.x, &scene.skyColor1.y, &scene.skyColor1.z,
            &scene.skyColor2.x, &scene.skyColor2.y, &scene.skyColor2.z);
    if (c != 7) {
        scene.gammaCorrection = 2.2f;
        vinit(scene.skyColor1, 0.15f, 0.3f, 0.5f);
        vinit(scene.skyColor2, 1.f, 1.f, 1.f);
    }

	// part 2: camera configuration: camera eye.x eye.y eye.z  target.x target.y target.z
	c = fscanf(f,"camera %f %f %f  %f %f %f\n",
			&camera.orig.x, &camera.orig.y, &camera.orig.z,
			&camera.target.x, &camera.target.y, &camera.target.z);
	camera.pitch = 0;
	camera.yaw = 0;
	camera.width = width;
	camera.height = height;

	originalCamera = camera;
	if (c != 6) {
		fprintf(stderr, "Failed to read 6 camera parameters: %d\n", c);
		exit(-1);
	}

	// part 3: material counts: materials n
	c = fscanf(f, "materials %u\n", &scene.materialCount);
	if (c < 1) {
        fprintf(stderr, "Failed to read the number of materials: %d\n", c);
        exit(-1);
	}
	// creates space for the 3 default materials
	#define NUMBER_OF_DEFAULT_MATERIALS 3
	scene.materialCount += NUMBER_OF_DEFAULT_MATERIALS;

	// part 4: material descriptors
	materials = (Material *)malloc(sizeof(Material) * scene.materialCount);
	materials[0].type = LAMBERTIAN;
	materials[1].type = CONDUCTOR;
	materials[2].type = DIELECTRIC;

	unsigned int lineIndex;
	for (lineIndex = NUMBER_OF_DEFAULT_MATERIALS; lineIndex < scene.materialCount; lineIndex++) {
        Material *material = &materials[lineIndex];
        char materialTypeString[100];

        int readValuesCount = fscanf(f, "%s", materialTypeString);
        if (readValuesCount != 1) {
            fprintf(stderr, "Failed to read material description type. Expected 1, found %d value(s)\n", readValuesCount);
            exit(-1);
        }

        if (strcmp(materialTypeString, "matte") == 0) {
            material->type = MATTE;
            readValuesCount = fscanf(f,"%f %f %f  %f\n",
                &material->kd.x, &material->kd.y, &material->kd.z,
                &material->sigma);
            if (readValuesCount != 4) {
                fprintf(stderr, "Failed to read matte material descriptor. Expected 4, found %d value(s)\n", readValuesCount);
                exit(-1);
            }
            clamp(material->sigma, 0.f, 0.9f);
//            material->sigma = ((float)FLOAT_PI/180.f) * (material->sigma);
            float sigma2 = material->sigma*material->sigma;
            material->A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
            material->B = 0.45f * sigma2 / (sigma2 + 0.09f);

        } else if (strcmp(materialTypeString, "plastic") == 0) {
            material->type = PLASTIC;
            readValuesCount = fscanf(f,"%f %f %f  %f %f %f  %f\n",
                &material->kd.x, &material->kd.y, &material->kd.z,
                &material->ks.x, &material->ks.y, &material->ks.z,
                &material->roughness);
            if (readValuesCount != 7) {
                fprintf(stderr, "Failed to read plastic material descriptor. Expected 4, found %d value(s)\n", readValuesCount);
                exit(-1);
            }
        }
	}

	// part 5: object counts: objects n
	c = fscanf(f,"objects %u\n", &objectCount);
	if (c < 1) {
		fprintf(stderr, "Failed to read object count: %d\n", c);
		exit(-1);
	}
	fprintf(stderr, "Scene size: %d\n", objectCount);

	// part 6: object descriptors
	objects = (Object *)malloc(sizeof(Object) * objectCount);
	unsigned int addedObjects = 0;
	lightCount = 0;
	for (lineIndex = 0; lineIndex < objectCount; lineIndex++) {
        Object *obj = &objects[lineIndex + addedObjects];
		unsigned int mat;
		char objectTypeString[100];

		int readValuesCount = fscanf(f, "%s", objectTypeString);
        if (readValuesCount != 1) {
            fprintf(stderr, "Failed to read object description type. Expected 1, found %d value(s)\n", readValuesCount);
            exit(-1);
        }
//        printf("reading an object of type %s and saving on position %d.\n", objectTypeString, lineIndex + addedObjects);

		if (strcmp(objectTypeString, "sphere") == 0) {
		    obj->type = SPHERE;
            // line x: sphere     radius  xyz  emission  color  reflection
            readValuesCount = fscanf(f,"%f  %f %f %f  %f %f %f  %f %f %f  %u\n",
				&obj->radius,
				&obj->center.x, &obj->center.y, &obj->center.z,
				&obj->emission.x, &obj->emission.y, &obj->emission.z,
				&obj->color.x, &obj->color.y, &obj->color.z,
				&mat);

            if (readValuesCount != 11) {
                fprintf(stderr, "Failed to read sphere #%d: %d\n", lineIndex, readValuesCount);
                exit(-1);
            }

            obj->materialId = getMaterialType(mat, lineIndex);
            printf("Sphere material id: %d\n", obj->materialId);
            obj->area = 4 * FLOAT_PI * obj->radius*obj->radius;

            if (!viszero(obj->emission)) {
                lightCount++;
            }
		} else if (strcmp(objectTypeString, "triangle") == 0) {
		    obj->type = TRIANGLE;
		    Object *tri = obj;
            readValuesCount = fscanf(f, "%f %f %f  %f %f %f  %f %f %f  %f %f %f  %f %f %f  %d\n",
                &obj->p1.x, &obj->p1.y, &obj->p1.z,
                &obj->p2.x, &obj->p2.y, &obj->p2.z,
                &obj->p3.x, &obj->p3.y, &obj->p3.z,
                &obj->emission.x, &obj->emission.y, &obj->emission.z,
                &obj->color.x, &obj->color.y, &obj->color.z,
                &mat);

            if (readValuesCount != 16) {
                fprintf(stderr, "Failed to read triangle #%d: %d\n", lineIndex, readValuesCount);
                exit(-1);
            }
//            printf("triangle: %f %f %f   %f %f %f   %f %f %f   %f %f %f   %f %f %f    %d\n",
//                   obj->p1.x, obj->p1.y, obj->p1.z, obj->p2.x, obj->p2.y, obj->p2.z, obj->p3.x, obj->p3.y, obj->p3.z, obj->emission.x, obj->emission.y, obj->emission.z, obj->color.x, obj->color.y, obj->color.z, obj->refl);
            obj->materialId = getMaterialType(mat, lineIndex);


            // the area of the circle outside the triangle is (abc)/4area
            float a, b, c;
            a = dist(tri->p2, tri->p1);
            b = dist(tri->p3, tri->p2);
            c = dist(tri->p1, tri->p3);
            vec e1; vsub(e1, tri->p2, tri->p1);
            vec e2; vsub(e2, tri->p3, tri->p1);
            vec normal; vxcross(normal, e1, e2);
            tri->area = norm(normal) * 0.5f;
            vnorm(normal);

            printf("tri #%d normal: %.2f %.2f %.2f\n", lineIndex + addedObjects, normal.x, normal.y, normal.z);


            if (!viszero(tri->emission)) {
                lightCount++;
//                printf("Tri emitindo luz: %f %f %f\t%f %f %f\t%f %f %f\n", tri->p1.x, tri->p1.y, tri->p1.z, tri->p2.x, tri->p2.y, tri->p2.z, tri->p3.x, tri->p3.y, tri->p3.z);
//                printf("Color: %f %f %f\tEmission: %f %f %f\n", tri->color.x, tri->color.y, tri->color.z, tri->emission.x, tri->emission.y, tri->emission.z);
//                printf("Area: %f\n", tri->area);
//                printf("Raio: %f\n", tri->radius);
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

            obj->materialId = getMaterialType(mat, lineIndex);

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

            unsigned int flags = (unsigned int) NULL;
            int ret = tinyobj_parse_obj(&attrib, &shapes, &num_shapes, &materials,
                                    &num_materials, data, data_len, flags);
            if (ret != TINYOBJ_SUCCESS) {
                fprintf(stderr, "Error loading obj. Exiting...");
                exit(1);
            }

//            printf("# of shapes    = %d\n", (int)num_shapes);
//            printf("# of materials = %d\n", (int)num_materials);


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
                tri->materialId = obj->materialId;

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

                // the area of the circle outside the triangle is (abc)/4area
                float a, b, c;
                a = dist(tri->p2, tri->p1);
                b = dist(tri->p3, tri->p2);
                c = dist(tri->p1, tri->p3);
                vec e1; vsub(e1, tri->p2, tri->p1);
                vec e2; vsub(e2, tri->p3, tri->p1);
                vec normal; vxcross(normal, e1, e2);
                tri->area = norm(normal) * 0.5f;


                if (!viszero(tri->emission)) {
                    lightCount++;
//                    printf("Tri emitindo luz: %f %f %f\t%f %f %f\t%f %f %f\n", p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
//                    printf("Color: %f %f %f\tEmission: %f %f %f\n", tri->color.x, tri->color.y, tri->color.z, tri->emission.x, tri->emission.y, tri->emission.z);
                }


                faceOffset += 3;
            }

            addedObjects += numTriangles;
		}

	}

	objectCount += addedObjects;

	printf("Finished parsing scene. Resulting objectCount = %d\n", objectCount);

//	int i;
//	for (int i = 0; i < objectCount; i++) {
//        printf("Object #%d: %d\n", i, objects[i].type);
//	}

	fclose(f);

//	int *counts = (int*) malloc(sizeof(int) * 2);
//	counts[0] = objectCount;
//	counts[1] = lightCount;
//	return counts;
}

void updateCamera(int delta) {
    bool cameraMoved = false;
    float speed =  delta * 0.0001f * (shiftPressed ? 10 : 1);
    vec movement = {0.f, 0.f, 0.f};

    // forward or backward movement
    if (keyStates['w'] || keyStates['W'] || keyStates['s'] || keyStates['S']) {
        const int frontOrBack = keyStates['w'] || keyStates['W'] ? 1 : -1;
        vsmul(movement, frontOrBack, camera.dir);

        cameraMoved = true;
    }

    // lateral movement
    if (keyStates['a'] || keyStates['A'] || keyStates['d'] || keyStates['D']) {
        int leftOrRight = keyStates['d'] || keyStates['D'] ? 1 : -1;

        vec direction = camera.x;
        vnorm(direction);
        vsmul(direction, leftOrRight, direction);

        vadd(movement, direction, movement);
        cameraMoved = true;
    }

    // if there was movement, we scale it by the max speed and integrate it
    // on the camera position & target
    if (cameraMoved) {
        vnorm(movement);
        vsmul(movement, speed, movement);

        vadd(camera.orig, camera.orig, movement);
        vadd(camera.target, camera.target, movement);

        reInitViewpointDependentBuffers(0);
    }
}

///
/// updates the opengl window by asking opencl to execute the kernel and then
/// asking glut to redraw the window
///
void update(void) {
    int timeNow = glutGet(GLUT_ELAPSED_TIME);
    int delta = timeNow - (timeBefore || 0);
    updateCamera(delta);
	updateRendering();

	timeBefore = timeNow;
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

	// renders the stats at the bottom
	glColor3f(1.f, 1.f, 1.f);
	glRasterPos2i(4, 23);
	renderString(GLUT_BITMAP_8_BY_13, captionLine1);
	glRasterPos2i(4, 10);
	renderString(GLUT_BITMAP_8_BY_13, captionLine2);

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


void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	    case 'P':
		case 'p': {
		    char fileName[200];
		    char timeSinceBeginning[20];
		    getHumanReadableTime(wallClockTime() - startRenderingTime, timeSinceBeginning);
		    sprintf(fileName, "%s-%i-%s.ppm", sceneTitle, currentSample, timeSinceBeginning);
			FILE *f = fopen(fileName, "w"); // Write image to PPM file.
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
		case ' ':
            // resets the camera to its initial position
            camera = originalCamera;
            reInitViewpointDependentBuffers(0);
            break;

		case 27: // ESC
			fprintf(stderr, "Done.\n");
			exit(0);
			break;

        case 'w':
        case 'W':
            setKeyStates('w', 'W', true);
            break;

        case 's':
        case 'S':
            setKeyStates('s', 'S', true);
            break;

        case 'a':
        case 'A':
            setKeyStates('a', 'A', true);
            break;

        case 'd':
        case 'D':
            setKeyStates('d', 'D', true);
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
		default:
			break;
	}
	keyStates[key] = true;
}

inline void setKeyStates(unsigned char lowercase, unsigned char uppercase, bool value) {
    keyStates[lowercase] = value;
    keyStates[uppercase] = value;
}

/// simply records the key that was released
void keyboardUp(unsigned char key, int x, int y) {
    switch (key) {
    case 'w':
    case 'W':
        setKeyStates('w', 'W', false);
        break;

    case 's':
    case 'S':
        setKeyStates('s', 'S', false);
        break;

    case 'a':
    case 'A':
        setKeyStates('a', 'A', false);
        break;

    case 'd':
    case 'D':
        setKeyStates('d', 'D', false);
        break;
    }
}


/// gets the current mouse position and compares to the last one to check if
/// we need to change the pitch/yaw of the camera
/// it only changes the camera if the mouse button is pressed
void motion(int x, int y) {
	int deltaX = mouseX - x;
	int deltaY = mouseY - y;

	if (deltaX != 0 || deltaY != 0) {

        // rotate the camera using pitch (nodding movement) and yaw (nonono movement)
		if (mouseButton== GLUT_LEFT_BUTTON) {
			camera.yaw += deltaX * 0.01;
			camera.yaw = camera.yaw - TWO_PI * floor(camera.yaw / TWO_PI);
			camera.pitch += -deltaY * 0.01;
			camera.pitch = clamp(camera.pitch, -PI_OVER_TWO, PI_OVER_TWO);
		}

		glutSetCursor(GLUT_CURSOR_CROSSHAIR);

		mouseX = x;
		mouseY = y;
		reInitViewpointDependentBuffers(0);
	} else {
        glutSetCursor(GLUT_CURSOR_INHERIT);
    }
}

/// simply records the mouse state
void mouse(int button, int state, int x, int y) {
	mouseButton = button;
	mouseX = x;
	mouseY = y;

	motion(x, y);
}

void special(int key, int x, int y) {
    if (key == GLUT_KEY_SHIFT_L || key == GLUT_KEY_SHIFT_R) {
        shiftPressed = true;
    }
}

void specialUp(int key, int x, int y) {
    if (key == GLUT_KEY_SHIFT_L || key == GLUT_KEY_SHIFT_R) {
        shiftPressed = false;
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

	for (int i = 0; i < 256; i++) {
        keyStates[i] = false;
	}

    printf("About to register callbacks...\n");
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutKeyboardUpFunc(keyboardUp);
    glutSpecialFunc(special);
    glutSpecialUpFunc(specialUp);
    glutDisplayFunc(displayScene);
	glutIdleFunc(update);

    printf("About to config viewport and projection...\n");
	glViewport(0, 0, width, height);
	glLoadIdentity();
	glOrtho(0.f, width - 1.f, 0.f, height - 1.f, -1.f, 1.f);
}
