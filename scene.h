#ifndef _SCENE_H_
#define _SCENE_H_

#include "vec.h"

typedef struct {
    float gammaCorrection;
    vec skyColor1, skyColor2;
    unsigned int materialCount;
} Scene;

#endif // _SCENE_H_
