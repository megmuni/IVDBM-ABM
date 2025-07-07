/*
 * timer.h
 *
 * File Contents: Contains declarations for functions for timing code sections
 *
 * Created on: May 26, 2015
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

#include <unistd.h>
#include <sys/time.h>

using namespace std;

/*
 * Description:	Function to start timer
 *
 * Return: Time at the start
 *
 * Parameters: void
 */
struct timeval tic();

/*
 * Description:	Function to stop timer
 *
 * Return: void
 *
 * Parameters: timeval  -- Time at the start
 */
long toc(struct timeval);

/*
 * Description:	Output time info
 *
 * Return: void
 *
 * Parameters: time        -- Elapsed time in millisecond
 *             stageName   -- The name of the stage to be displayed
 *             stageLevel  -- The stage level to be displayed
 */
void print_time(
		long time,
		const char* stageName,
		const char* stageLevel);

/*
 * Description:	A macro for timing a command/function of a major stage in go() and print the timing info
 *
 * Return: void
 *
 * Parameters: task    -- Command/function to perform
 *             sName   -- The name of the stage to be displayed
 *             sLevel  -- The stage level to be displayed
 *
 * Example Usage:
 * 			TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");
 */
#define TIME_STAGE(task, sName, sLevel) do { \
	struct timeval START = tic ();\
	task;\
	long ET = toc(START);\
	print_time(ET, sName, sLevel);\
	} while (0)

#endif /* TIMER_H_ */
