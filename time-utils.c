#include "time-utils.h"
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

double wallClockTime() {
	struct timeval time;
	gettimeofday(&time, NULL);

	return time.tv_sec + time.tv_usec / 1000000.0;
}


/// returns (on the time param) a human readable string for a certain
/// number of seconds, eg, 15ms, 56s, 1m, 24h
/// numbers are always integers, rounded up or down
void getHumanReadableTime(const double seconds, char *time) {
    if (seconds < 1) {
        // less than 1s: show in ms
        sprintf(time, "%03.0fms", roundf(seconds * 1000.f));
    } else if (seconds < 1 * 60) {
        // less than 1min: show in s
        sprintf(time, "%.0fs", roundf(seconds));
    } else if (seconds < 1 * 60 * 60) {
        // less than 1h: show in m
        sprintf(time, "%.0fm%.0fs", floorf(seconds / 60.f), fmod(floorf(seconds), 60.f));
    } else {
        // 1h or more: show in h
        sprintf(time, "%.0fh%02.0fm", floorf(seconds / 60.f / 60.f), fmod(floorf(seconds / 60.f), 60.f));
    }
}

