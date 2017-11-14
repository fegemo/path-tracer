#ifndef _DISPLAYFUNC_H
#define	_DISPLAYFUNC_H

#include <math.h>
#include <GL/freeglut.h>

#include "vec.h"

extern int width;
extern int height;
extern unsigned int *pixels;
extern unsigned int renderingFlags;
extern char captionBuffer[256];

extern int amiSmallptCPU;

extern void initGlut(int argc, char *argv[], char *windowTittle);
extern double wallClockTime();

extern int readScene(char *);
extern void updateCamera();

#endif	/* _DISPLAYFUNC_H */

