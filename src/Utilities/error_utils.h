/*
 * error_utils.h
 * 
 * File Contents: Contains error reporting functions that can be used for 
 *                functions that return an error code
 *
 * Created on: May 27, 2015
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef ERROR_UTILS_H_
#define ERROR_UTILS_H_

#pragma once

#include <stdio.h>

namespace util {

/*
 * Description:	Displays error message if error is non-zero
 *
 * Return: The error message from the function called
 *
 * Parameters: error     -- Error code
 *             message   -- Error message
 *             filename	 -- Error source code file name
 *             line		   -- Error source line number
 *             print		 -- true if want to print the error (default)
 *
 * Example Usage:
 * if (ABMerror(f(param), "WoundHealingWorld f() failed", __FILE__, __LINE__ ))
 *  exit(1);
 */
int ABMerror(int error,
		const char* message,
		const char* filename,
		int line,
		bool print = true) {
	if (error && print) {
		fprintf(stderr, "[%s, %d] %s (WH ABMs error %d)\n", filename, line, message, error);
		fflush(stderr);
	}
	return error;
}

/*
 * Description:	Displays error message if error is non-zero
 *
 * Return: The error message from the function called
 *
 * Parameters: error  -- Error code
 *             print  -- true if want to print the error (default)
 */
int ABMerror(
		int error,
		bool print = true) {
	if (error && print) {
		fprintf(stderr, "(WH ABMs error %d)\n", error);
		fflush(stderr);
	}
	return error;
}

} //namespace util

#endif /* ERROR_UTILS_H_ */
