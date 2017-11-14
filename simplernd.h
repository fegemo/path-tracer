#ifndef _SIMPLERND_H
#define	_SIMPLERND_H

/*
 * A Simple Random number generator
 * from http://en.wikipedia.org/wiki/Random_number_generator
 */

#ifndef SMALLPT_GPU

static float GetRandom(unsigned int *seed0, unsigned int *seed1) {
	*seed0 = 36969 * ((*seed0) & 65535) + ((*seed0) >> 16);
	*seed1 = 18000 * ((*seed1) & 65535) + ((*seed1) >> 16);

	unsigned int ires = ((*seed0) << 16) + (*seed1);

	/* Convert to float */
	union {
		float f;
		unsigned int ui;
	} res;
	res.ui = (ires & 0x007fffff) | 0x40000000;

	return (res.f - 2.f) / 2.f;
}

#endif

#endif	/* _SIMPLERND_H */

