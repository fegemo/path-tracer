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

#ifndef _GEOMFUNC_H
#define	_GEOMFUNC_H

#include "geom.h"
#include "simplernd.h"

#ifndef SMALLPT_GPU

static float SphereIntersect(
    OCL_CONSTANT_BUFFER const Sphere *s,
	const Ray *r) { /* returns distance, 0 if nohit */
	vec op; /* Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 */
	vsub(op, s->p, r->o);

	float b = vdot(op, r->d);
	float det = b * b - vdot(op, op) + s->rad * s->rad;
	if (det < 0.f)
		return 0.f;
	else
		det = sqrt(det);

	float t = b - det;
	if (t >  EPSILON)
		return t;
	else {
		t = b + det;

		if (t >  EPSILON)
			return t;
		else
			return 0.f;
	}
}

///
/// finds a random point on a sphere
///  u1, u2 are random [0,1]
///  v is the output point
///
/// this is from pbrt page 664
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
/// checks if a ray intersects an object in the scene and return the intersection information
///
static int Intersect(
    OCL_CONSTANT_BUFFER const Sphere *spheres,
	const unsigned int sphereCount,
	const Ray *r,
	float *t,
	unsigned int *id) {

	float inf = (*t) = 1e20f;

	unsigned int i = sphereCount;
	for (; i--;) {
        // ASSUMES SPHERES
		const float d = SphereIntersect(&spheres[i], r);
		if ((d != 0.f) && (d < *t)) {
			*t = d;
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
    OCL_CONSTANT_BUFFER const Sphere *spheres,
	const unsigned int sphereCount,
	// the ray being cast
	const Ray *r,
	// the maximum distance
	const float maxt) {

	unsigned int i = sphereCount;

	// iterates on the scene checking for intersection with the ray
	for (; i--;) {
        // ASSUMES SPHERES
		const float d = SphereIntersect(&spheres[i], r);
		if ((d != 0.f) && (d < maxt))
			return 1;
	}

	return 0;
}

///
/// iterates over the lights in the scene to get how much each one contribute directly to
/// the point (hitPoint) from a surface with a (normal) and output to (result)
///
static void SampleLights(
	// the scene
	OCL_CONSTANT_BUFFER const Sphere *spheres,
	// the size of the scene
	const unsigned int sphereCount,
	// random seeds
	unsigned int *seed0, unsigned int *seed1,
	// the point we're calculating the direct contribution of 01 light
	const vec *hitPoint,
	// the normal on that point
	const vec *normal,
	// the returning direct lighting contribution
	vec *result) {

	// result = (0,0,0)
	vclr(*result);

	// for each light...
	unsigned int i;
	for (i = 0; i < sphereCount; i++) {
		OCL_CONSTANT_BUFFER const Sphere *light = &spheres[i];
		if (!viszero(light->e)) {
			// this is a light source (it has an emission component)...
			// the shadow ray starts at the hitPoint and goes to a random point on the light
			Ray shadowRay;
			shadowRay.o = *hitPoint;

			// choose a random point over the area of the light source
			// ASSUMES SPHERE (...i think... maybe this could be applicable to "bounding spheres"...)
			// first we find a random point in a unit sphere, then we multiply it by the light sphere radius
			// and then we find the point in world coordinates (by adding the light sphere center coords)
			vec unitSpherePoint;
			UniformSampleSphere(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &unitSpherePoint);
			vec spherePoint;
			vsmul(spherePoint, light->rad, unitSpherePoint);
			vadd(spherePoint, spherePoint, light->p);

			// configures the direction of the shadow ray, finds its length (keep it on len)
			// and then normalizes the direction vector
			vsub(shadowRay.d, spherePoint, *hitPoint);
			const float len = sqrt(vdot(shadowRay.d, shadowRay.d));
			vsmul(shadowRay.d, 1.f / len, shadowRay.d);

			// finds the angle between shadow and the light sphere normal
			// if it is accute (ie, dot > 0, we hit the backside of the light
            //    ...unitSpherePoint is a vector that points from the sphere center to its surface point
			float wo = vdot(shadowRay.d, unitSpherePoint);
			if (wo > 0.f) {
				// we hit the other half of the sphere... should ignore it
				continue;
			} else {
			    // we just flip the sign of the cosine because we want it positive
			    //   ...basically because we want the vector to be from the light sphere surface to its center
				wo = -wo;
			}

			// now we send the shadow ray to the light and see if it hits another object before it...
			//   wi > 0 <=> the light is in the outside and in the direction of the hitPoint normal
			const float wi = vdot(shadowRay.d, *normal);
			if ((wi > 0.f) && (!IntersectP(spheres, sphereCount, &shadowRay, len - EPSILON))) {
                // the object faces the light (wi > 0) and the shadow ray doesnt intersect anything
                //
                //  s = 4PIr² * cos(shadowR and light incidence) * cos(shadowR and normal) / d²
                //     ...where did this formula come from?
				vec lightColor; vassign(lightColor, light->e);
				const float s = (4.f * FLOAT_PI * light->rad * light->rad) * wi * wo / (len *len);
				vsmul(lightColor, s, lightColor);

				// last, we add the direct contribution of this light to the output Ld radiance
				vadd(*result, *result, lightColor);
			}
		}
	}
}

static void RadiancePathTracing(
    // the scene
    OCL_CONSTANT_BUFFER const Sphere *spheres,
	// size of the scene
	const unsigned int sphereCount,
	// primary ray
	const Ray *startRay,
	// seeds for random numbers
	unsigned int *seed0, unsigned int *seed1,
	// the radiance found for the ray
	vec *result) {

	Ray currentRay; rassign(currentRay, *startRay);
	vec radiance; vinit(radiance, 0.f, 0.f, 0.f);
	vec throughput; vinit(throughput, 1.f, 1.f, 1.f);

	unsigned int depth = 0;
	int specularBounce = 1;
	for (;; ++depth) {
		// Removed Russian Roulette in order to improve execution on SIMT
		// we bounce the ray for 6 times (primary + 6 = 7)
		if (depth > 6) {
			*result = radiance;
			return;
		}

		// distance to intersection
		float t;
		// id of intersected object (ie, its index in array)
		unsigned int id = 0;

		if (!Intersect(spheres, sphereCount, &currentRay, &t, &id)) {
			*result = radiance;
			// if the ray missed, return just the color so far...
			// we could put a skybox here...
			return;
		}

		// the hit object
		OCL_CONSTANT_BUFFER const Sphere *obj = &spheres[id];

		// the intersection point
		vec hitPoint;
		vsmul(hitPoint, t, currentRay.d);
		vadd(hitPoint, currentRay.o, hitPoint);

		// the normal at the intersection point
		// ASSUMING SPHERE
		vec normal;
		vsub(normal, hitPoint, obj->p);
		vnorm(normal);

		// cosine of the angle between the normal and the ray direction
		// this angle is accute if hitting from inside the sphere, and obtuse otherwise
		const float dp = vdot(normal, currentRay.d);

		// nl is the normal of the hitPoint and here we check
		// if we need to flip it (in case we're hitting from inside the sphere)
		// ASSUMING SPHERE
		vec nl;
		const float invSignDP = -1.f * sign(dp);
		vsmul(nl, invSignDP, normal);

		// adds the light emitted by the object hit
		vec eCol; vassign(eCol, obj->e);
		// if the object does emit light....
		if (!viszero(eCol)) {
			// and if this is a bounce originated from specular materials...
			if (specularBounce) {
				// reduce the emitted color according to the ray hitting angle and normal, then by this throughput,
				// then adds it to the output radiance
				vsmul(eCol, fabs(dp), eCol);
				vmul(eCol, throughput, eCol);
				vadd(radiance, radiance, eCol);
			}

			// light emitting objects do not bounce the ray, so we return
			*result = radiance;
			return;
		}

		// 100% DIFFUSE material
		if (obj->refl == DIFF) {
			// should not do the specular bounce... reduces the throughput by the object's color (whatever this is)
			specularBounce = 0;
			vmul(throughput, throughput, obj->c);

			// direct lighting component (splitting)
			// -------------------------------------

			vec Ld;
			SampleLights(spheres, sphereCount, seed0, seed1, &hitPoint, &nl, &Ld);
			vmul(Ld, throughput, Ld);
			vadd(radiance, radiance, Ld);

			// diffuse component (would be recursive, made iterative for perf.)
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
			vmul(throughput, throughput, obj->c);

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
				vmul(throughput, throughput, obj->c);

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
				vmul(throughput, throughput, obj->c);

				rassign(currentRay, reflRay);
				continue;
			} else {
				vsmul(throughput, TP, throughput);
				vmul(throughput, throughput, obj->c);

				rinit(currentRay, hitPoint, transDir);
				continue;
			}
		}
	}
}

static void RadianceDirectLighting(
	OCL_CONSTANT_BUFFER const Sphere *spheres,
	const unsigned int sphereCount,
	const Ray *startRay,
	unsigned int *seed0, unsigned int *seed1,
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

		float t; /* distance to intersection */
		unsigned int id = 0; /* id of intersected object */
		if (!Intersect(spheres, sphereCount, &currentRay, &t, &id)) {
			*result = rad; /* if miss, return */
			return;
		}

		OCL_CONSTANT_BUFFER const Sphere *obj = &spheres[id]; /* the hit object */

		vec hitPoint;
		vsmul(hitPoint, t, currentRay.d);
		vadd(hitPoint, currentRay.o, hitPoint);

		vec normal;
		vsub(normal, hitPoint, obj->p);
		vnorm(normal);

		const float dp = vdot(normal, currentRay.d);

		vec nl;
		// SIMT optimization
		const float invSignDP = -1.f * sign(dp);
		vsmul(nl, invSignDP, normal);

		/* Add emitted light */
		vec eCol; vassign(eCol, obj->e);
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
			vmul(throughput, throughput, obj->c);

			/* Direct lighting component */

			vec Ld;
			SampleLights(spheres, sphereCount, seed0, seed1, &hitPoint, &nl, &Ld);
			vmul(Ld, throughput, Ld);
			vadd(rad, rad, Ld);

			*result = rad;
			return;
		} else if (obj->refl == SPEC) { /* Ideal SPECULAR reflection */
			specularBounce = 1;

			vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			vmul(throughput, throughput, obj->c);

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
				vmul(throughput, throughput, obj->c);

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
				vmul(throughput, throughput, obj->c);

				rassign(currentRay, reflRay);
				continue;
			} else {
				vsmul(throughput, TP, throughput);
				vmul(throughput, throughput, obj->c);

				rinit(currentRay, hitPoint, transDir);
				continue;
			}
		}
	}
}

#endif

#endif	/* _GEOMFUNC_H */

