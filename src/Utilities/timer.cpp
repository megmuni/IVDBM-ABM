/*
 * timer.cpp
 *
 * File Contents: Contains functions to time sections of code 
 *
 * Created on: May 26, 2015
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#include "timer.h"

struct timeval tic () {
	struct timeval start;
	gettimeofday(&start, NULL);
	return start;
}

long toc (struct timeval start) {
	long out;
	static struct timeval end;
	gettimeofday(&end, NULL);
	out = (end.tv_sec*1000 + end.tv_usec/1000) - (start.tv_sec*1000 + start.tv_usec/1000);
	return out;  // in milliseconds
}

void print_time(long time, const char* stageName, const char* stageLevel) {
	printf("\tStage %s:\t%s took %ld ms\n", stageLevel, stageName, time);
}
