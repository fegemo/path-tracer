#ifndef _GEOMFUNC_H
#define	_GEOMFUNC_H

#include "geom.h"
#include "simplernd.h"

#ifndef SMALLPT_GPU

static float SphereIntersect(
    OCL_GLOBAL_BUFFER const Object *s,
	const Ray *r) { /* returns distance, 0 if nohit */
	vec op; /* Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 */
	vsub(op, s->center, r->o);

	float b = vdot(op, r->d);
	float delta = b * b - vdot(op, op) + s->radius * s->radius;
	if (delta < 0.f) {
		return 0.f;
	}
	else {
		delta = sqrt(delta);
	}

	float t = b - delta;
	if (t >  EPSILON) {
		return t;
	}
	else {
		t = b + delta;

		if (t >  EPSILON) {
			return t;
		}
		else {
			return 0.f;
		}
	}
}

///
/// Determines if a ray intersected with a triangle (algorithm from pbrt: 141-2)
/// Returns the t != 0 if there was a hit.
///
static float TriangleIntersect(OCL_GLOBAL_BUFFER const Object *triangle, const Ray *ray) {

    // calculates s1 = ray.direction x e2
    vec e1; vsub(e1, triangle->p2, triangle->p1);
    vec e2; vsub(e2, triangle->p3, triangle->p1);
    vec s1; vxcross(s1, ray->d, e2);
    float divisor = vdot(s1, e1);
    if (divisor == 0) {
        return 0;
    }
    float invDivisor = 1.f / divisor;

    // finds the first barycentric coordinate (b1)
    vec d; vsub(d, ray->o, triangle->p1);
    float b1 = vdot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.) {
        return 0;
    }

    // finds the second one (b2)
    vec s2; vxcross(s2, d, e1);
    float b2 = vdot(ray->d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.) {
        return 0;
    }

    // computes t to intersection point
    float t = vdot(e2, s2) * invDivisor;
    if (t > EPSILON) {
        return t;
    }

    return 0;
}

static void SphereNormal(vec *normal, OCL_GLOBAL_BUFFER const Object *sphere, const vec hitPoint) {
    vsub(*normal, hitPoint, sphere->center);
    vnorm(*normal);
}

static void TriangleNormal(vec *normal, OCL_GLOBAL_BUFFER const Object *triangle, vec hitPoint) {
    vec e1; vsub(e1, triangle->p2, triangle->p1);
    vec e2; vsub(e2, triangle->p3, triangle->p1);
    vxcross(*normal, e1, e2);
}

///
/// finds a random point on a sphere
///  u1, u2 are random [0,1]
///  v is the output point
///
/// this is from pbrt page 664
///
static void UniformSampleSphere(const float u1, const float u2, vec *v) {
    // zz = [-1, 1]
	const float zz = 1.f - 2.f * u1;
	// r = [0, 1] (trending towards 1)
	const float r = sqrt(max(0.f, 1.f - zz * zz));
	// phi = [0, 2pi]
	const float phi = 2.f * FLOAT_PI * u2;
	// xx =
	const float xx = r * cos(phi);
	const float yy = r * sin(phi);

	vinit(*v, xx, yy, zz);
}

///
/// finds a random point on a triangle
///  u1, u2 are random [0,1]
///  (u,v) is the barycentric coordinates for the output point
///
static void UniformSampleTriangle(const float u1, const float u2, float *u, float *v) {
    float su1 = (float)sqrt(u1);
    *u = 1.0f - su1;
    *v = u2 * su1;
}

///
/// checks if a ray intersects an object in the scene and returns the intersection information
///
static int Intersect(
    OCL_GLOBAL_BUFFER const Object *objects,
	const unsigned int objectCount,
	const Ray *ray,
	float *t,
	unsigned int *id) {

	float inf = (*t) = 1e20f;

	for (unsigned int i = 0; i < objectCount; i++) {
        float distance;
        switch (objects[i].type) {
        case SPHERE:
            distance = SphereIntersect(&objects[i], ray);
            break;
        case TRIANGLE:
            distance = TriangleIntersect(&objects[i], ray);
            break;
        case MODEL:
        default:
            continue;
        }

		if ((distance != 0.f) && (distance < *t)) {
			*t = distance;
			*id = i;
		}
	}


	return (*t < inf);
}

///
/// checks if a ray intersects the the scene without calculating hitpoint, normal etc.
///
static int IntersectP(
    // scene description
    OCL_GLOBAL_BUFFER const Object *objects,
	const unsigned int objectCount,
	// the ray being cast
	const Ray *ray,
	// the maximum distance
	const float maxt) {

	// iterates on the scene checking for intersection with the ray
	for (unsigned int i = 0; i < objectCount; i++) {
		float distance;

		switch (objects[i].type) {
        case SPHERE:
            distance = SphereIntersect(&objects[i], ray);
            break;
        case TRIANGLE:
            distance = TriangleIntersect(&objects[i], ray);
            break;
        case MODEL:
        default:
            continue;
		}

		if ((distance != 0.f) && (distance < maxt)) {
			return 1;
		}
	}

	return 0;
}

///
/// iterates over the lights in the scene to get how much each one contribute directly to
/// the point (hitPoint) from a surface with a (normal) and output to (result)
///
static void SampleLights(
	// the scene
	OCL_GLOBAL_BUFFER const Object *objects,
	// the size of the scene
	const unsigned int objectCount,
	const unsigned int lightCount,
	// random seeds
	unsigned int *seed0, unsigned int *seed1,
	// the point we're calculating the direct contribution from light
	const vec *hitPoint,
	// the normal on that point
	const vec *normal,
	// the returning direct lighting contribution
	vec *result) {

	// result = (0,0,0)
	vclr(*result);

	// for each light...
	unsigned int i;
	unsigned int lightsVisited;
	for (i = 0, lightsVisited = 0; i < objectCount && lightsVisited < lightCount; i++) {
		OCL_GLOBAL_BUFFER const Object *light = &objects[i];
		if (light->type == MODEL) continue;
		if (!viszero(light->emission)) {
			// this is a light source (as it has an emission component)...
            lightsVisited++;

			// the shadow ray starts at the hitPoint and goes to a random point on the light
			Ray shadowRay;
			shadowRay.o = *hitPoint;

			// choose a random point over the area of the light source
			vec lightPoint;
			vec lightNormalWhereShadowRayHit;
			switch (light->type) {
            case TRIANGLE:
                {
                    float barycentricU, barycentricV;
                    vec e1, e2;
                    UniformSampleTriangle(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &barycentricU, &barycentricV);
                    vsub(e1, light->p2, light->p1)
                    vsub(e2, light->p3, light->p1)
                    vsmul(e1, barycentricU, e1);
                    vsmul(e2, barycentricV, e2);
                    vassign(lightPoint, light->p1);
                    vadd(lightPoint, lightPoint, e1);
                    vadd(lightPoint, lightPoint, e2);

                    // sets what is the normal on the triangle where the shadow ray hit it
                    // (it's the same regardless of the point because it's a triangle)
                    vxcross(lightNormalWhereShadowRayHit, e1, e2);
                    vnorm(lightNormalWhereShadowRayHit);
                }
                break;

            case SPHERE:
            default:
                {
                    vec unitSpherePoint;
                    UniformSampleSphere(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &unitSpherePoint);
                    vsmul(lightPoint, light->radius, unitSpherePoint);
                    vadd(lightPoint, lightPoint, light->center);

                    // sets what is the normal on the triangle where the shadow ray hit it
                    vassign(lightNormalWhereShadowRayHit, unitSpherePoint);
                }
			}

			// configures the direction of the shadow ray, finds its length (keep it on len)
			// and then normalizes the direction vector
			vsub(shadowRay.d, lightPoint, *hitPoint);
			const float len = sqrt(vdot(shadowRay.d, shadowRay.d));
			vsmul(shadowRay.d, 1.f / len, shadowRay.d);

			// finds the angle between shadow and the light normal
			// if it is accute (ie, dot > 0, we hit the backside of the light
			float wo = vdot(shadowRay.d, lightNormalWhereShadowRayHit);
			if (wo > 0.f) {
				// we hit the other half of the light... should ignore it
				continue;
			} else {
			    // we just flip the sign of the cosine because we want it positive
			    //   ...basically because we want the vector to be from the light sphere surface to its center
				wo = -wo;
			}
			wo = clamp(wo, 0.f, 1.f);

			// now we send the shadow ray to the light and see if it hits another object before it...
			//   wi > 0 <=> the light is in the outside and in the direction of the hitPoint normal
			const float wi = vdot(shadowRay.d, *normal);
			if ((wi > 0.f) && (!IntersectP(objects, objectCount, &shadowRay, len - EPSILON))) {
                // the object faces the light (wi > 0) and the shadow ray doesnt intersect anything
                //
				vec lightColor; vassign(lightColor, light->emission);
                //  s = 4PIr² * cos(shadowR and light incidence) * cos(shadowR and normal) / d²
				const float s = light->area * wi * wo / (len *len);
				vsmul(lightColor, s, lightColor);

				// last, we add the direct contribution of this light to the output Ld radiance
				vadd(*result, *result, lightColor);
			}
		}
	}
	// if I take the average of the lights contributions the scene gets too dark...
	// although that's how pbrt does it...
    //vsmul(*result, 1.f/lightCount, *result);
}

///
/// randomly chooses one light in the scene to get how much it directly contributes to
/// the point (hitPoint) from a surface with a (normal) and output to (result)
///
//static void UniformSampleOneLight(
//	// the scene
//	OCL_GLOBAL_BUFFER const Object *objects,
//	// the size of the scene
//	const unsigned int objectCount,
//	const unsigned int lightCount,
//	// random seeds
//	unsigned int *seed0, unsigned int *seed1,
//	// the point we're calculating the direct contribution from light
//	const vec *hitPoint,
//	// the normal on that point
//	const vec *normal,
//	// the returning direct lighting contribution
//	vec *result) {
//
//	// result = (0,0,0)
//	vclr(*result);
//
//	// chooses one light
//	unsigned int lightIndex = (int)floor(GetRandom(seed0, seed1) * lightCount);
//	// for each light...
//	unsigned int i;
//	unsigned int lightsVisited;
//	for (i = 0, lightsVisited = 0; i < objectCount && lightsVisited < lightCount; i++) {
//		OCL_GLOBAL_BUFFER const Object *light = &objects[i];
//		if (light->type == MODEL) continue;
//		if (viszero(light->emission)) continue;
//        if (lightsVisited++ != lightIndex) continue;
//
//        // the shadow ray starts at the hitPoint and goes to a random point on the light
//        Ray shadowRay;
//        shadowRay.o = *hitPoint;
//
//        // choose a random point over the area of the light source
//        vec lightPoint;
//        vec lightNormalWhereShadowRayHit;
//        switch (light->type) {
//        case TRIANGLE:
//            {
//                float barycentricU, barycentricV;
//                vec e1, e2;
//                UniformSampleTriangle(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &barycentricU, &barycentricV);
//                vsub(e1, light->p2, light->p1)
//                vsub(e2, light->p3, light->p1)
//                vsmul(e1, barycentricU, e1);
//                vsmul(e2, barycentricV, e2);
//                vassign(lightPoint, light->p1);
//                vadd(lightPoint, lightPoint, e1);
//                vadd(lightPoint, lightPoint, e2);
//
//                // sets what is the normal on the triangle where the shadow ray hit it
//                // (it's the same regardless of the point because it's a triangle)
//                vxcross(lightNormalWhereShadowRayHit, e1, e2);
//                vnorm(lightNormalWhereShadowRayHit);
//            }
//            break;
//
//        case SPHERE:
//        default:
//            {
//
//                vec unitSpherePoint;
//                UniformSampleSphere(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &unitSpherePoint);
//                vsmul(lightPoint, light->radius, unitSpherePoint);
//                vadd(lightPoint, lightPoint, light->center);
//
//                // sets what is the normal on the triangle where the shadow ray hit it
//                vassign(lightNormalWhereShadowRayHit, unitSpherePoint);
//            }
//        }
//
//        // configures the direction of the shadow ray, finds its length (keep it on len)
//        // and then normalizes the direction vector
//        vsub(shadowRay.d, lightPoint, *hitPoint);
//        const float len = sqrt(vdot(shadowRay.d, shadowRay.d));
//        vsmul(shadowRay.d, 1.f / len, shadowRay.d);
//
//        // finds the angle between shadow and the light normal
//        // if it is accute (ie, dot > 0, we hit the backside of the light
//        float wo = vdot(shadowRay.d, lightNormalWhereShadowRayHit);
//        if (wo > 0.f) {
//            // we hit the other half of the light... should ignore it
//            continue;
//        } else {
//            // we just flip the sign of the cosine because we want it positive
//            //   ...basically because we want the vector to be from the light sphere surface to its center
//            wo = -wo;
//        }
//        wo = clamp(wo, 0.f, 1.f);
//
//        // now we send the shadow ray to the light and see if it hits another object before it...
//        //   wi > 0 <=> the light is in the outside and in the direction of the hitPoint normal
//        const float wi = vdot(shadowRay.d, *normal);
//        if ((wi > 0.f) && (!IntersectP(objects, objectCount, &shadowRay, len - EPSILON))) {
//            // the object faces the light (wi > 0) and the shadow ray doesnt intersect anything
//            //
//            vec lightColor; vassign(lightColor, light->emission);
//            //  s = 4PIr² * cos(shadowR and light incidence) * cos(shadowR and normal) / d²
//            const float s = light->area * wi * wo / (len *len);
//            vsmul(lightColor, s, lightColor);
//
//            // last, we add the direct contribution of this light to the output Ld radiance
//            vadd(*result, *result, lightColor);
//        }
//
//        // exits the for, as we only wanted to sample from 01 light and we did it
//        vsmul(*result, lightCount, *result);
//        break;
//	}
//}



static vec ComputeBackgroundColor(const vec direction, const vec sky1, const vec sky2) {
    vec skyDirection = {-.41f, -.82f, .41f};
    float position = (vdot(direction, skyDirection) + 1) / 2;
	vec interpolatedColor; affine(interpolatedColor, position, sky1, sky2);
	return interpolatedColor;
}

static void RadiancePathTracing(
    // the scene
    OCL_GLOBAL_BUFFER const Object *objects,
	// size of the scene (objects and lights)
	const unsigned int objectCount,
	const unsigned int lightCount,
	// primary ray
	const Ray *startRay,
	// seeds for random numbers
	unsigned int *seed0, unsigned int *seed1,
	// sky colors
	const vec sky1, const vec sky2,
	// the radiance found for the ray
	vec *result) {

	Ray currentRay; rassign(currentRay, *startRay);
	vec rad; vinit(rad, 0.f, 0.f, 0.f);
	vec throughput; vinit(throughput, 1.f, 1.f, 1.f);

	unsigned int depth = 0;
	int specularBounce = 1;
	for (;; ++depth) {
		// no russian roulette so threads finish at the same time
		// we bounce the ray for 6 times (primary + 6 = 7)
		if (depth > 6) {
			*result = rad;
			return;
		}

		// distance to intersection
		float t;
		// id of intersected object (ie, its index in array)
		unsigned int id = 0;

		if (!Intersect(objects, objectCount, &currentRay, &t, &id)) {
			// if the ray missed, return a sky color
            vec backgroundColor = ComputeBackgroundColor(currentRay.d, sky1, sky2);
            vmul(*result, throughput, backgroundColor);

			return;
		}

		// the hit object
		OCL_GLOBAL_BUFFER const Object *obj = &objects[id];

		// the intersection point
		vec hitPoint;
		vsmul(hitPoint, t, currentRay.d);
		vadd(hitPoint, currentRay.o, hitPoint);

		// the normal at the intersection point
		vec normal;
		switch (obj->type) {
        case SPHERE:
            SphereNormal(&normal, obj, hitPoint);
            break;

        case TRIANGLE:
            TriangleNormal(&normal, obj, hitPoint);
            break;
		}
		vnorm(normal);

		// cosine of the angle between the normal and the ray direction
		// this angle is accute if hitting from inside the sphere, and obtuse otherwise
		const float dp = vdot(normal, currentRay.d);

		// nl is the normal of the hitPoint and here we check
		// if we need to flip it (in case we're hitting from inside the sphere)
		vec nl;
		const float invSignDP = -1.f * sign(dp);
		vsmul(nl, invSignDP, normal);

		// adds the light emitted by the object hit
		vec emittedColor; vassign(emittedColor, obj->emission);
		// if the object does emit light....
		if (!viszero(emittedColor)) {
			// and if this is a bounce originated from specular materials...
			if (specularBounce) {
				// reduce the emitted color according to the ray hitting angle and normal, then by this throughput,
				// then adds it to the output radiance
				vsmul(emittedColor, fabs(dp), emittedColor);
				vmul(emittedColor, throughput, emittedColor);
				vadd(rad, rad, emittedColor);
			}

			// light emitting objects do not bounce the ray, so we return
			*result = rad;
			return;
		}

		// 100% DIFFUSE material
		if (obj->refl == DIFF) { /* Ideal DIFFUSE reflection */
			// should not do the specular bounce... reduces the throughput by the object's color (whatever this is)
			specularBounce = 0;
			vmul(throughput, throughput, obj->color);

			// direct lighting component (splitting)
			// -------------------------------------

			vec Ld;
			SampleLights(objects, objectCount, lightCount, seed0, seed1, &hitPoint, &nl, &Ld);
			//UniformSampleOneLight(objects, objectCount, lightCount, seed0, seed1, &hitPoint, &nl, &Ld);
			vmul(Ld, throughput, Ld);
			vadd(rad, rad, Ld);

			// diffuse component
			// ----------------------------------------------------------------

			// we'll sample a position in a unit circle that will be used to find how this currentRay
			// will be reflected (as this is a diffuse material, it can go uniformly in any direction
            // on the hemisphere)
			float r1 = 2.f * FLOAT_PI * GetRandom(seed0, seed1);
			float r2 = GetRandom(seed0, seed1);
			float r2s = sqrt(r2);

			// w = the normal of the hitpoint
			vec w; vassign(w, nl);

			// u, v, w = orthonormal basis on the hitpoint of the object ('a' is just temporary)
			vec u, a;
			if (fabs(w.x) > .1f) {
				vinit(a, 0.f, 1.f, 0.f);
			} else {
				vinit(a, 1.f, 0.f, 0.f);
			}
			vxcross(u, a, w);
			vnorm(u);

			vec v;
			vxcross(v, w, u);

			// finds the direction of the next ray on our path to trace
			vec newDir;
			vsmul(u, cos(r1) * r2s, u);
			vsmul(v, sin(r1) * r2s, v);
			vadd(newDir, u, v);
			vsmul(w, sqrt(1 - r2), w);
			vadd(newDir, newDir, w);

			currentRay.o = hitPoint;
			currentRay.d = newDir;

			// finished this ray, let's go shoot the the next one
			continue;
		} else if (obj->refl == SPEC) {
			// 100% SPECULAR material
			specularBounce = 1;

			// finds the perfectly reflected ray
			vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			// multiplies the throughput by the object color (whatever that is...)
			vmul(throughput, throughput, obj->color);

			// sets up the next ray in the series and shoots it in the scene
			rinit(currentRay, hitPoint, newDir);
			continue;
		} else {
			// 100% REFRACTION material
			specularBounce = 1;

			vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			Ray reflRay; rinit(reflRay, hitPoint, newDir); /* Ideal dielectric REFRACTION */
			int into = (vdot(normal, nl) > 0); /* Ray from outside going in? */

			float nc = 1.f;
			float nt = 1.5f;
			float nnt = into ? nc / nt : nt / nc;
			float ddn = vdot(currentRay.d, nl);
			float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

			if (cos2t < 0.f)  { /* Total internal reflection */
				vmul(throughput, throughput, obj->color);

				rassign(currentRay, reflRay);
				continue;
			}

			float kk = (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t));
			vec nkk;
			vsmul(nkk, kk, normal);
			vec transDir;
			vsmul(transDir, nnt, currentRay.d);
			vsub(transDir, transDir, nkk);
			vnorm(transDir);

			float a = nt - nc;
			float b = nt + nc;
			float R0 = a * a / (b * b);
			float c = 1 - (into ? -ddn : vdot(transDir, normal));

			float Re = R0 + (1 - R0) * c * c * c * c*c;
			float Tr = 1.f - Re;
			float P = .25f + .5f * Re;
			float RP = Re / P;
			float TP = Tr / (1.f - P);

			if (GetRandom(seed0, seed1) < P) { /* R.R. */
				vsmul(throughput, RP, throughput);
				vmul(throughput, throughput, obj->color);

				rassign(currentRay, reflRay);
				continue;
			} else {
				vsmul(throughput, TP, throughput);
				vmul(throughput, throughput, obj->color);

				rinit(currentRay, hitPoint, transDir);
				continue;
			}
		}
	}
}

static void RadianceDirectLighting(
	OCL_GLOBAL_BUFFER const Object *objects,
	const unsigned int objectCount,
	const unsigned int lightCount,
	const Ray *startRay,
	unsigned int *seed0, unsigned int *seed1,
	// sky colors
	const vec sky1, const vec sky2,
	vec *result) {

	Ray currentRay; rassign(currentRay, *startRay);
	vec rad; vinit(rad, 0.f, 0.f, 0.f);
	vec throughput; vinit(throughput, 1.f, 1.f, 1.f);

	unsigned int depth = 0;
	int specularBounce = 1;

	for (;; ++depth) {
		// Removed Russian Roulette in order to improve execution on SIMT
		if (depth > 6) {
			*result = rad;
			return;
		}

		float t;
		unsigned int id = 0;
		if (!Intersect(objects, objectCount, &currentRay, &t, &id)) {
			// if the ray missed, return a sky color
            vec backgroundColor = ComputeBackgroundColor(currentRay.d, sky1, sky2);
            vmul(*result, throughput, backgroundColor);

			return;
		}

		OCL_GLOBAL_BUFFER const Object *obj = &objects[id]; /* the hit object */

		vec hitPoint;
		vsmul(hitPoint, t, currentRay.d);
		vadd(hitPoint, currentRay.o, hitPoint);

		// the normal at the intersection point
		vec normal;
		switch (obj->type) {
        case SPHERE:
            SphereNormal(&normal, obj, hitPoint);
            break;

        case TRIANGLE:
            TriangleNormal(&normal, obj, hitPoint);
            break;

		}
		vnorm(normal);

		const float dp = vdot(normal, currentRay.d);

		vec nl;
		// SIMT optimization
		const float invSignDP = -1.f * sign(dp);
		vsmul(nl, invSignDP, normal);

		/* Add emitted light */
		vec eCol; vassign(eCol, obj->emission);
		if (!viszero(eCol)) {
			if (specularBounce) {
				vsmul(eCol, fabs(dp), eCol);
				vmul(eCol, throughput, eCol);
				vadd(rad, rad, eCol);
			}

			*result = rad;
			return;
		}

		if (obj->refl == DIFF) { /* Ideal DIFFUSE reflection */
			specularBounce = 0;
			vmul(throughput, throughput, obj->color);

			/* Direct lighting component */

			vec Ld;
			SampleLights(objects, objectCount, lightCount, seed0, seed1, &hitPoint, &nl, &Ld);
			vmul(Ld, throughput, Ld);
			vadd(rad, rad, Ld);

			*result = rad;
			return;
		} else if (obj->refl == SPEC) { /* Ideal SPECULAR reflection */
			specularBounce = 1;

			vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			vmul(throughput, throughput, obj->color);

			rinit(currentRay, hitPoint, newDir);
			continue;
		} else {
			specularBounce = 1;

			vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			Ray reflRay; rinit(reflRay, hitPoint, newDir); /* Ideal dielectric REFRACTION */
			int into = (vdot(normal, nl) > 0); /* Ray from outside going in? */

			float nc = 1.f;
			float nt = 1.5f;
			float nnt = into ? nc / nt : nt / nc;
			float ddn = vdot(currentRay.d, nl);
			float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

			if (cos2t < 0.f)  { /* Total internal reflection */
				vmul(throughput, throughput, obj->color);

				rassign(currentRay, reflRay);
				continue;
			}

			float kk = (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t));
			vec nkk;
			vsmul(nkk, kk, normal);
			vec transDir;
			vsmul(transDir, nnt, currentRay.d);
			vsub(transDir, transDir, nkk);
			vnorm(transDir);

			float a = nt - nc;
			float b = nt + nc;
			float R0 = a * a / (b * b);
			float c = 1 - (into ? -ddn : vdot(transDir, normal));

			float Re = R0 + (1 - R0) * c * c * c * c*c;
			float Tr = 1.f - Re;
			float P = .25f + .5f * Re;
			float RP = Re / P;
			float TP = Tr / (1.f - P);

			if (GetRandom(seed0, seed1) < P) { /* R.R. */
				vsmul(throughput, RP, throughput);
				vmul(throughput, throughput, obj->color);

				rassign(currentRay, reflRay);
				continue;
			} else {
				vsmul(throughput, TP, throughput);
				vmul(throughput, throughput, obj->color);

				rinit(currentRay, hitPoint, transDir);
				continue;
			}
		}
	}
}

#endif

#endif	/* _GEOMFUNC_H */

