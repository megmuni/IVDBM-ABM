/*
 * woundHealingWorld.cpp
 *
 * File contents: Contains the WHWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <tgmath.h>
#define PI 3.14159
#include <string>
#include <sstream>
#include <omp.h>
#include "woundHealingWorld.h"
#include "../../enums.h"
#include "../../Utilities/timer.h"
#include "../../Utilities/error_utils.h"
#include "../../Utilities/input_utils.h"

//Helper functions for CUDA
#include <helper_functions.h>
#include <helper_cuda.h>
#include <helper_timer.h>

/*
	// Include CUDA runtime and CUFFT
	#include <cuda_runtime.h>
	#include <cufft.h>
	//Helper functions for CUDA
	#include <helper_functions.h>
	#include <helper_cuda.h>
	#include "../../Diffusion/convolutionFFT2D_common.h"
	#include "../../Diffusion/test.cuh"
*/

#define KCoeffsH		3
#define KCoeffsW		3

typedef struct CCTX		// convolution context
{
	int KH;
	int KW;
	int KX;
	int KY;
	int DH;
	int DW;
	int FFTH;
	int FFTW;
	int windowH;
	int windowW;
} c_ctx;

using namespace std;

#ifdef GPU_DIFFUSE	// (*)

	/* -------------------------------------------------------------------------- */
	/*                              Helper functions                              */
	/* -------------------------------------------------------------------------- */
	int snapTransformSize(int dataSize)
	{
		int hiBit;
		unsigned int lowPOT, hiPOT;
		
		dataSize = iAlignUp(dataSize, 16);
		for (hiBit = 31; hiBit >= 0; hiBit--) if (dataSize & (1U << hiBit)) break;
			
		lowPOT = 1U << hiBit;
		if (lowPOT == (unsigned int)dataSize) return dataSize;

		hiPOT = 1U << (hiBit + 1);
		if (hiPOT <= 1024) return hiPOT;
		else return pow(2, ceil(log(dataSize)/log(2))); //return iAlignUp(dataSize, 512);
	}

	bool computeKernel_Old(
			float		*h_ResultGPU,
			float		*d_UnpaddedResult,
			int			 kernelRadius,
			float		 lambda,
			float		 gamma,					// decay constant
			float		 dt,
			c_ctx 		 cctx)
	{

		float t = 0.0;

		int
		cpu_input  = 0,
		cpu_output = 1;

		float
		*h_Data,
		*h_Kernel,
		*h_ResultCPU[2];

		float
		*d_Data,
		*d_PaddedData,
		*d_Kernel,
		*d_PaddedKernel;

		fComplex
		*d_DataSpectrum,
		*d_KernelSpectrum;

		cufftHandle
		fftPlanFwd,
		fftPlanInv;

		bool bRetVal;
		StopWatchInterface *hTimer = NULL;
		sdkCreateTimer(&hTimer);

		fprintf(stderr, "Testing kernel computation\n");
		printf("Testing kernel computation\n");
		fprintf(stderr, "\tBuilding filter kernel\n");
		const int    kernelH = cctx.KH;//7;
		const int    kernelW = cctx.KW;//6;
		const int    kernelY = cctx.KY;//3;
		const int    kernelX = cctx.KX;//4;
		const int      dataH = cctx.DH;//100;//1160;//2000;
		const int      dataW = cctx.DW;//100;//1660;//2000;
		const int outKernelH = cctx.DH;
		const int outKernelW = cctx.DW;
		const int       fftH = cctx.FFTH;
		const int       fftW = cctx.FFTW;

		fprintf(stderr, "...allocating memory\n");
		h_Data      			= (float *)malloc(dataH   * dataW * sizeof(float));
		h_Kernel    			= (float *)malloc(kernelH * kernelW * sizeof(float));
		h_ResultCPU[cpu_input] 	= (float *)malloc(dataH   * dataW * sizeof(float));
		h_ResultCPU[cpu_output] = (float *)malloc(dataH   * dataW * sizeof(float));

		checkCudaErrors(cudaMalloc((void **)&d_Data,   dataH   * dataW   * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_Kernel, kernelH * kernelW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum,   fftH * (fftW / 2 + 1) * sizeof(fComplex)));
		checkCudaErrors(cudaMalloc((void **)&d_KernelSpectrum, fftH * (fftW / 2 + 1) * sizeof(fComplex)));
		fprintf(stderr, "...generating 2D %d x %d kernel coefficients\n", kernelH, kernelW);

		for (int i = 0; i < kernelH * kernelW; i++) h_Kernel[i] = 0;
		h_Kernel[0 * kernelW + 1] = lambda;
		h_Kernel[1 * kernelW + 0] = lambda;
		h_Kernel[1 * kernelW + 2] = lambda;
		h_Kernel[2 * kernelW + 1] = lambda;
		h_Kernel[1 * kernelW + 1] = 1 - 4*lambda - gamma*dt;

		for (int i = 0; i < dataH * dataW; i++){
			h_Data[i] = 0;
			h_ResultCPU[cpu_input][i] = 0;
		}

		// Copy kernel data to middle block of the input
		int start_i = outKernelH/2 - kernelH/2;
		int end_i	= outKernelH/2 + kernelH/2 + 1;
		int start_j = outKernelW/2 - kernelW/2;
		int end_j	= outKernelW/2 + kernelW/2 + 1;
		int ki = 0, kj = 0;
		for (int i = start_i; i < end_i; i++) {
			for (int j = start_j; j < end_j; j++) {
				h_Data					[i * dataW + j]	= h_Kernel[ki * kernelW + kj];
				h_ResultCPU [cpu_input]	[i * dataW + j]	= h_Kernel[ki * kernelW + kj];
				printf("%d,%d -> %d,%d\n", ki, kj, i, j);
				kj++;
			}
			ki++;
			kj = 0;
		}

		fprintf(stderr, "...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
		checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
		checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));
		fprintf(stderr, "...uploading to GPU and padding convolution kernel and input data\n");
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaMemcpy(d_Kernel, h_Kernel, kernelH * kernelW * sizeof(float), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_Data,   h_Data,   dataH   * dataW *   sizeof(float), cudaMemcpyHostToDevice));
		sdkStopTimer(&hTimer);
		
		double dataTransferTime = sdkGetTimerValue(&hTimer);
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMemset(d_PaddedData,   0, fftH * fftW * sizeof(float)));

		padKernel(
				d_PaddedKernel,
				d_Kernel,
				fftH,
				fftW,
				kernelH,
				kernelW,
				kernelY,
				kernelX
		);

		padDataClampToBorder(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX
		);

		sdkStopTimer(&hTimer);
		double memsetPaddingTime = sdkGetTimerValue(&hTimer);

		//Not including kernel transformation into time measurement, since convolution kernel is not changed very frequently
		fprintf(stderr, "...transforming convolution kernel\n");

		double buildKernelTimeTotalGPU = 0;
		double buildKernelTimeTotalCPU = 0;
		checkCudaErrors(cudaDeviceSynchronize());

		// d_KernelSpectrum = FFT(d_PaddedKernel)
		checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));
		for (int iter = 0; iter < kernelRadius; iter++){
			fprintf(stderr, "...running GPU Kernel building iteration %d:\n", iter);
			printf("GPU Kernel building iteration %d:\n", iter);

			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);
			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));
			modulateAndNormalize(d_DataSpectrum, d_KernelSpectrum, fftH, fftW, 1);
			checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));
			checkCudaErrors(cudaDeviceSynchronize());
			sdkStopTimer(&hTimer);
			
			double gpuTime = sdkGetTimerValue(&hTimer);
			printf("\t\tGPU computation: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (gpuTime * 0.001), gpuTime);
			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);
			fprintf(stderr, "...removing results padding\n");
		/*	
			unpadResult(
					d_UnpaddedResult,
					d_PaddedData,
					dataH,
					dataW,
					fftH,
					fftW);
		*/
			sdkStopTimer(&hTimer);
			double unpadTime = sdkGetTimerValue(&hTimer);
			printf("\t\textract results: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (unpadTime * 0.001), unpadTime);
			fprintf(stderr, "...reading back GPU convolution results\n");
			checkCudaErrors(cudaMemcpy(h_ResultGPU, d_UnpaddedResult, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));
			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);
			fprintf(stderr, "...running reference CPU convolution\n");
	//		convolutionClampToBorderCPU(
	//				h_ResultCPU[cpu_output],
	//				h_ResultCPU[cpu_input],
	//				h_Kernel,
	//				dataH,
	//				dataW,
	//				kernelH,
	//				kernelW,
	//				kernelY,
	//				kernelX
	//		);

			sdkStopTimer(&hTimer);
			double cpuTime = sdkGetTimerValue(&hTimer);
			printf("\t\tCPU computation: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (cpuTime * 0.001), cpuTime);

			buildKernelTimeTotalGPU += gpuTime;
			buildKernelTimeTotalGPU += unpadTime;
			buildKernelTimeTotalCPU += cpuTime;

			t += dt;
			// Update indices for CPU input/output
			cpu_input  = (cpu_input  + 1) % 2;
			cpu_output = (cpu_output + 1) % 2;
		}

		for (int i = 0; i < kernelH; i++) {
			for (int j = 0; j < kernelW; j++) {
				printf(", %f", h_Kernel[i*kernelW + j]);
			}
			printf("\n");
		}

		printf("\tData transfer:		%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
		printf("\tMemset and padding:	%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (memsetPaddingTime * 0.001), memsetPaddingTime);
		printf("\tTotal GPU time:		%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalGPU * 0.001), buildKernelTimeTotalGPU);
		printf("\tTotal CPU time:		%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalCPU * 0.001), buildKernelTimeTotalCPU);
		fprintf(stderr, "...comparing the results: ");
		double sum_delta2 = 0;
		double sum_ref2   = 0;
		double max_delta_ref = 0;

		// Update indices for CPU input/output
		cpu_input  = (cpu_input  + 1) % 2;
		cpu_output = (cpu_output + 1) % 2;

		fprintf(stderr, "...shutting down\n");
		sdkDeleteTimer(&hTimer);
		checkCudaErrors(cufftDestroy(fftPlanInv));
		checkCudaErrors(cufftDestroy(fftPlanFwd));
		checkCudaErrors(cudaFree(d_DataSpectrum));
		checkCudaErrors(cudaFree(d_KernelSpectrum));
		checkCudaErrors(cudaFree(d_PaddedKernel));
		checkCudaErrors(cudaFree(d_Data));
		checkCudaErrors(cudaFree(d_Kernel));
		checkCudaErrors(cudaFree(d_PaddedData));

		free(h_ResultCPU[0]);
		free(h_ResultCPU[1]);
		free(h_Data);
		free (h_Kernel);

		return bRetVal;
	}

	bool computeKernelSpectrum_Old(
		fComplex	*d_KernelSpectrum,
		float		*d_Kernel,
		c_ctx 		kernel_cctx,
		c_ctx		chem_cctx
		){

		float
		*d_UnpaddedKernel,
		*d_PaddedKernel;

		cufftHandle
		fftPlanFwd,
		fftPlanInv;

		bool bRetVal;

		printf("Testing kernel spectrum computation\n");
		fprintf(stderr, "Testing kernel spectrum computation\n");
		const int    kernelH = chem_cctx.KH;//kernel_cctx.DH;//7;
		const int    kernelW = chem_cctx.KW;//kernel_cctx.DW;//6;
		const int    kernelY = chem_cctx.KY;//kernel_cctx.DH / 2;//3;
		const int    kernelX = chem_cctx.KX;//kernel_cctx.DW / 2;//4;
		const int       fftH = chem_cctx.FFTH;
		const int       fftW = chem_cctx.FFTW;

		printf("\tkernelH: %d\tkernelW: %d\n", kernelH, kernelW);
		printf("\tkernelX: %d\tkernelY: %d\n", kernelX, kernelY);
		printf("\tfftH: %d\tfftW: %d\n", fftH, fftW);
		fprintf(stderr,"...allocating memory\n");
		checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_UnpaddedKernel, kernelH * kernelW * sizeof(float)));
		fprintf(stderr,"...creating R2C FFT plans for %i x %i\n", fftH, fftW);
		checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
		checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));
		fprintf(stderr,"...uploading to GPU and padding convolution kernel and input data\n");
		checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fftH * fftW * sizeof(float)));

		padKernel(
			d_PaddedKernel,
			d_Kernel,
			fftH,
			fftW,
			kernelH,
			kernelW,
			kernelY,
			kernelX
		);

		checkCudaErrors(cudaDeviceSynchronize());
		// d_KernelSpectrum = FFT(d_PaddedKernel)
		checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel,       (cufftComplex *)d_KernelSpectrum));
		checkCudaErrors(cudaDeviceSynchronize());
		checkCudaErrors(cufftDestroy(fftPlanFwd));
		checkCudaErrors(cudaFree(d_PaddedKernel));
		checkCudaErrors(cudaFree(d_UnpaddedKernel));
		return bRetVal;
	}

	bool computeChemDiffusionGPU_Old(
		float		*h_ChemOut,
		float		*h_ChemIn,
		fComplex	*d_KernelSpectrum,
		c_ctx		 cctx,
		int		 	 iter){

		float t = 0.0;

		float
		*d_Data,
		*d_PaddedData,
		*d_UnpaddedResult;

		fComplex
		*d_DataSpectrum;

		cufftHandle
		fftPlanFwd,
		fftPlanInv;

		bool bRetVal;
		StopWatchInterface *hTimer = NULL;
		sdkCreateTimer(&hTimer);

		printf("Testing GPU chemical diffusion computation\n");
		fprintf(stderr,"Testing GPU chemical diffusion computation\n");
		const int    kernelH = cctx.KH;
		const int    kernelW = cctx.KW;
		const int    kernelY = cctx.KY;
		const int    kernelX = cctx.KX;
		const int      dataH = cctx.DH;
		const int      dataW = cctx.DW;
		const int       fftH = cctx.FFTH;
		const int       fftW = cctx.FFTW;

		printf("\tkernelH: %d\tkernelW: %d\n", kernelH, kernelW);
		printf("\tkernelX: %d\tkernelY: %d\n", kernelX, kernelY);
		printf("\tdataH: %d\tdataW: %d\n", dataH, dataH);
		printf("\tfftH: %d\tfftW: %d\n", fftH, fftW);
		fprintf(stderr,"...allocating memory\n");
		checkCudaErrors(cudaMalloc((void **)&d_Data,   dataH   * dataW   * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_UnpaddedResult, dataH * dataW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum,   fftH * (fftW / 2 + 1) * sizeof(fComplex)));
		fprintf(stderr,"...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
		checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
		checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));
		fprintf(stderr,"...uploading to GPU and padding input data\n");
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaMemcpy(d_Data,   h_ChemIn,   dataH   * dataW *   sizeof(float), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemset(d_PaddedData,   0, fftH * fftW * sizeof(float)));
		sdkStopTimer(&hTimer);
		double dataTransferTime = sdkGetTimerValue(&hTimer);

	/*	padDataRightWall(
			d_PaddedData,
			d_Data,
			fftH,
			fftW,
			dataH,
			dataW,
			kernelH,
			kernelW,
			kernelY,
			kernelX
		);
	*/
		fprintf(stderr,"...performing convolution\n");
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaDeviceSynchronize());

		// --------- Computing convolution ------------ begin

		// d_DataSpectrum = FFT(d_PaddedData)
		checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));

		// d_DataSpectrum = d_DataSpectrum * d_KernelSpectrum
		modulateAndNormalize(d_DataSpectrum, d_KernelSpectrum, fftH, fftW, 1);

		// d_PaddedData = IFFT(d_DataSpectrum)			<------- Output
		checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));

		// --------- Computing convolution ------------ end

		checkCudaErrors(cudaDeviceSynchronize());
		sdkStopTimer(&hTimer);
		
		double gpuTime = sdkGetTimerValue(&hTimer);
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		fprintf(stderr,"...removing results padding\n");
	/*	unpadResult(
			d_UnpaddedResult,
			d_PaddedData,
			dataH,
			dataW,
			fftH,
			fftW
		);
	*/
		sdkStopTimer(&hTimer);
		double unpadTime = sdkGetTimerValue(&hTimer);
		fprintf(stderr,"...reading back GPU convolution results\n");
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		// h_ChemOut = d_UnpaddedResult
		checkCudaErrors(cudaMemcpy(h_ChemOut, d_UnpaddedResult, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));
		sdkStopTimer(&hTimer);
		double readbackTime = sdkGetTimerValue(&hTimer);

		printf("\t\tdata transfer:				%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
		printf("\t\tGPU chemical diffusion computation:	%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (gpuTime * 0.001), gpuTime);
		printf("\t\textract results:			%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (unpadTime * 0.001), unpadTime);
		printf("\t\tread back:				%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (readbackTime * 0.001), readbackTime);
		
		fprintf(stderr,"...shutting down\n");
		sdkDeleteTimer(&hTimer);
		checkCudaErrors(cufftDestroy(fftPlanInv));
		checkCudaErrors(cufftDestroy(fftPlanFwd));

		printf("...freeing device pointers\n");
		checkCudaErrors(cudaFree(d_DataSpectrum));
		checkCudaErrors(cudaFree(d_UnpaddedResult));
		checkCudaErrors(cudaFree(d_PaddedData));
		checkCudaErrors(cudaFree(d_Data));
		printf("...returning to main()\n");
		return bRetVal;
	}

	/*
		bool computeChemDiffusionCPU(
				float		*h_ChemOut,
				float		*h_ChemIn,
				float		*h_Kernel,
				c_ctx		 cctx,
				int		 	 iter){

			bool bRetVal = 1;
			StopWatchInterface *hTimer = NULL;
			sdkCreateTimer(&hTimer);
			printf("Testing Chemical Diffusion CPU\n");
			const int    kernelH = cctx.KH;//7;
			const int    kernelW = cctx.KW;//6;
			const int    kernelY = cctx.KY;//3;
			const int    kernelX = cctx.KX;//4;
			const int      dataH = cctx.DH;//100;//1160;//2000;
			const int      dataW = cctx.DW;//100;//1660;//2000;
			const int outKernelH = cctx.DH;
			const int outKernelW = cctx.DW;
			const int       fftH = cctx.FFTH;
			const int       fftW = cctx.FFTW;
			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);
			fprintf(stderr,"...running reference CPU convolution\n");
			convolutionClampToBorderCPU(
					h_ChemOut,
					h_ChemIn,
					h_Kernel,
					dataH,
					dataW,
					kernelH,
					kernelW,
					kernelY,
					kernelX
			);
			sdkStopTimer(&hTimer);
			double cpuTime = sdkGetTimerValue(&hTimer);
			printf("\t\tCPU chemical diffusion computation:\t%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (cpuTime * 0.001), cpuTime);
			sdkDeleteTimer(&hTimer);
			return bRetVal;
		}
	*/

	bool computeKernel(
			float		*h_ResultGPU,
			float		*d_UnpaddedResult,
			// Changed 2
			float		*h_Window,
			float		*d_Window,
			int			 kernelRadius,
			float		 lambda,
			float		 gamma,					// decay constant
			float		 dt,
			c_ctx 		 cctx){

		float t = 0.0;

		int
		cpu_input  = 0,
		cpu_output = 1;

		float
		*h_Data,
		*h_Kernel,
		*h_ResultCPU[2];

		float
		*d_Data,
		*d_PaddedData,
		*d_Kernel,
		*d_PaddedKernel;

		fComplex
		*d_DataSpectrum,
		*d_KernelSpectrum;

		cufftHandle
		fftPlanFwd,
		fftPlanInv;

		bool bRetVal;
		StopWatchInterface *hTimer = NULL;
		sdkCreateTimer(&hTimer);

		#ifdef PRINT_KERNEL
			fprintf(stderr, "Testing kernel computation\n");
			printf("Testing kernel computation\n");
			fprintf(stderr, "\tBuilding filter kernel\n");
		#endif

		const int    kernelH = cctx.KH;//7;
		const int    kernelW = cctx.KW;//6;
		const int    kernelY = cctx.KY;//3;
		const int    kernelX = cctx.KX;//4;
		const int      dataH = cctx.DH;//100;//1160;//2000;
		const int      dataW = cctx.DW;//100;//1660;//2000;
		const int outKernelH = cctx.DH;
		const int outKernelW = cctx.DW;
		const int       fftH = cctx.FFTH;
		const int       fftW = cctx.FFTW;
		// Changed 2
		const int    windowH = cctx.windowH;
		const int    windowW = cctx.windowW;

		#ifdef PRINT_KERNEL
			fprintf(stderr, "...allocating memory\n");
		#endif
		h_Data      			= (float *)malloc(dataH   * dataW * sizeof(float));
		h_Kernel    			= (float *)malloc(kernelH * kernelW * sizeof(float));
		h_ResultCPU[cpu_input] 	= (float *)malloc(dataH   * dataW * sizeof(float));
		h_ResultCPU[cpu_output] = (float *)malloc(dataH   * dataW * sizeof(float));
		checkCudaErrors(cudaMalloc((void **)&d_Data,   dataH   * dataW   * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_Kernel, kernelH * kernelW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum,   fftH * (fftW / 2 + 1) * sizeof(fComplex)));
		checkCudaErrors(cudaMalloc((void **)&d_KernelSpectrum, fftH * (fftW / 2 + 1) * sizeof(fComplex)));

		#ifdef PRINT_KERNEL
			fprintf(stderr, "...generating 2D %d x %d kernel coefficients\n", kernelH, kernelW);
		#endif

		for (int i = 0; i < kernelH * kernelW; i++) h_Kernel[i] = 0;
		
		h_Kernel[0 * kernelW + 1] = lambda;
		h_Kernel[1 * kernelW + 0] = lambda;
		h_Kernel[1 * kernelW + 2] = lambda;
		h_Kernel[2 * kernelW + 1] = lambda;
		h_Kernel[1 * kernelW + 1] = 1 - 4*lambda - gamma*dt;

		for (int i = 0; i < dataH * dataW; i++){
			h_Data[i] = 0;
			h_ResultCPU[cpu_input][i] = 0;
		}

		// Copy kernel data to middle block of the input
		int start_i = outKernelH/2 - kernelH/2;
		int end_i	= outKernelH/2 + kernelH/2 + 1;
		int start_j = outKernelW/2 - kernelW/2;
		int end_j	= outKernelW/2 + kernelW/2 + 1;
		int ki = 0, kj = 0;
		for (int i = start_i; i < end_i; i++) {
			for (int j = start_j; j < end_j; j++) {
				h_Data					[i * dataW + j]	= h_Kernel[ki * kernelW + kj];
				h_ResultCPU [cpu_input]	[i * dataW + j]	= h_Kernel[ki * kernelW + kj];
				#ifdef PRINT_KERNEL
					printf("%d,%d -> %d,%d\n", ki, kj, i, j);
				#endif
				kj++;
			}
			ki++;
			kj = 0;
		}

		#ifdef PRINT_KERNEL
			fprintf(stderr, "...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
		#endif
		checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
		checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));

		#ifdef PRINT_KERNEL
			fprintf(stderr, "...uploading to GPU and padding convolution kernel and input data\n");
		#endif
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaMemcpy(d_Kernel, h_Kernel, kernelH * kernelW * sizeof(float), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_Data,   h_Data,   dataH   * dataW *   sizeof(float), cudaMemcpyHostToDevice));
		sdkStopTimer(&hTimer);
		
		double dataTransferTime = sdkGetTimerValue(&hTimer);
		#ifdef PRINT_KERNEL
		//	printf("\tData transfer: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
		#endif
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMemset(d_PaddedData,   0, fftH * fftW * sizeof(float)));

		padKernel(
				d_PaddedKernel,
				d_Kernel,
				fftH,
				fftW,
				kernelH,
				kernelW,
				kernelY,
				kernelX
		);

		padDataClampToBorder(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX
		);

		sdkStopTimer(&hTimer);
		double memsetPaddingTime = sdkGetTimerValue(&hTimer);
		#ifdef PRINT_KERNEL
		//	printf("\tMemset and padding: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (memsetPaddingTime * 0.001), memsetPaddingTime);
		#endif

		//Not including kernel transformation into time measurement, since convolution kernel is not changed very frequently
		#ifdef PRINT_KERNEL
			fprintf(stderr, "...transforming convolution kernel\n");
		#endif

		double buildKernelTimeTotalGPU = 0;
		double buildKernelTimeTotalCPU = 0;
		checkCudaErrors(cudaDeviceSynchronize());
		// d_KernelSpectrum = FFT(d_PaddedKernel)
		checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));
		
		for (int iter = 0; iter < kernelRadius; iter++){

			#ifdef PRINT_KERNEL
				fprintf(stderr, "...running GPU Kernel building iteration %d:\n", iter);
				printf("GPU Kernel building iteration %d:\n", iter);
			#endif
			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);
			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));
			modulateAndNormalize(d_DataSpectrum, d_KernelSpectrum, fftH, fftW, 1);
			checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));
			checkCudaErrors(cudaDeviceSynchronize());
			sdkStopTimer(&hTimer);
			
			double gpuTime = sdkGetTimerValue(&hTimer);
			#ifdef PRINT_KERNEL
					printf("\t\tGPU computation: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (gpuTime * 0.001), gpuTime);
			#endif

			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);
			#ifdef PRINT_KERNEL
				fprintf(stderr, "...removing results padding\n");
			#endif
	/*		unpadResult(
					d_UnpaddedResult,
					d_PaddedData,
					dataH,
					dataW,
					fftH,
					fftW
			);
	*/
			sdkStopTimer(&hTimer);
			double unpadTime = sdkGetTimerValue(&hTimer);
	#ifdef PRINT_KERNEL
			printf("\t\textract results: %f MPix/s (%f ms)\n",
					(double)dataH * (double)dataW * 1e-6 / (unpadTime * 0.001), unpadTime);

			fprintf(stderr, "...reading back GPU convolution results\n");
	#endif
			checkCudaErrors(cudaMemcpy(h_ResultGPU, d_UnpaddedResult, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));

			sdkResetTimer(&hTimer);
			sdkStartTimer(&hTimer);

	#ifdef PRINT_KERNEL
			fprintf(stderr, "...running reference CPU convolution\n");
	#endif

	// convolutionClampToBorderCPU(
	//				h_ResultCPU[cpu_output],
	//				h_ResultCPU[cpu_input],
	//				h_Kernel,
	//				dataH,
	//				dataW,
	//				kernelH,
	//				kernelW,
	//				kernelY,
	//				kernelX
	//		);

			sdkStopTimer(&hTimer);
			double cpuTime = sdkGetTimerValue(&hTimer);
			#ifdef PRINT_KERNEL
				printf("\t\tCPU computation: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (cpuTime * 0.001), cpuTime);
			#endif

			buildKernelTimeTotalGPU += gpuTime;
			buildKernelTimeTotalGPU += unpadTime;
			buildKernelTimeTotalCPU += cpuTime;
			t += dt;

			// Update indices for CPU input/output
			cpu_input  = (cpu_input  + 1) % 2;
			cpu_output = (cpu_output + 1) % 2;
		}

		// Changed 2 : Added
		#ifdef PRINT_KERNEL
			printf("...extract kernel window from center\n");
		#endif
	/*	extractCenter(
				d_Window,
				d_UnpaddedResult,
				dataH,
				dataW,
				windowH,
				windowW
		);
	*/

		// Changed 2 : Added
		#ifdef PRINT_KERNEL
			printf("...reading back kernel center from GPU\n");
		#endif
		checkCudaErrors(cudaMemcpy(h_Window, d_Window, windowH * windowW * sizeof(float), cudaMemcpyDeviceToHost));

		for (int i = 0; i < kernelH; i++) {
			for (int j = 0; j < kernelW; j++) {
				#ifdef PRINT_KERNEL
					printf(", %f", h_Kernel[i*kernelW + j]);
				#endif
			}
			#ifdef PRINT_KERNEL
				printf("\n");
			#endif
		}

		#ifdef PRINT_KERNEL
			printf("\tData transfer:		%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
			printf("\tMemset and padding:	%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (memsetPaddingTime * 0.001), memsetPaddingTime);
			printf("\tTotal GPU time:		%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalGPU * 0.001), buildKernelTimeTotalGPU);
			printf("\tTotal CPU time:		%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalCPU * 0.001), buildKernelTimeTotalCPU);
			fprintf(stderr, "...comparing the results: ");
		#endif
		double sum_delta2 = 0;
		double sum_ref2   = 0;
		double max_delta_ref = 0;

		// Update indices for CPU input/output
		cpu_input  = (cpu_input  + 1) % 2;
		cpu_output = (cpu_output + 1) % 2;

	//	for (int y = 0; y < dataH; y++)
	//		for (int x = 0; x < dataW; x++){
	//			double  rCPU = (double)h_ResultCPU[cpu_output][y * dataW + x];
	//			double  rGPU = (double)h_ResultGPU[y * dataW  + x];
	//			double delta = (rCPU - rGPU) * (rCPU - rGPU);
	//			double   ref = rCPU * rCPU + rCPU * rCPU;
	//
	//			if ((delta / ref) > max_delta_ref) max_delta_ref = delta / ref;
	//			sum_delta2 += delta;
	//			sum_ref2   += ref;
	//		}
	//
	//	double L2norm = sqrt(sum_delta2 / sum_ref2);
	//	printf("rel L2 = %E (max delta = %E)\n", L2norm, sqrt(max_delta_ref));
	//	bRetVal = (L2norm < 1e-6) ? true : false;
	//	printf(bRetVal ? "L2norm Error OK\n" : "L2norm Error too high!\n");

	#ifdef PRINT_KERNEL
		fprintf(stderr, "...shutting down\n");
	#endif
		sdkDeleteTimer(&hTimer);
		checkCudaErrors(cufftDestroy(fftPlanInv));
		checkCudaErrors(cufftDestroy(fftPlanFwd));
		checkCudaErrors(cudaFree(d_DataSpectrum));
		checkCudaErrors(cudaFree(d_KernelSpectrum));
		checkCudaErrors(cudaFree(d_PaddedKernel));
		checkCudaErrors(cudaFree(d_Data));
		checkCudaErrors(cudaFree(d_Kernel));
		checkCudaErrors(cudaFree(d_PaddedData));
		free(h_ResultCPU[0]);
		free(h_ResultCPU[1]);
		free(h_Data);
		free(h_Kernel);

		return bRetVal;
	}


	bool computeKernelSpectrum(
		fComplex	*d_KernelSpectrum,
		float		*d_Kernel,
		c_ctx 		kernel_cctx,
		c_ctx		chem_cctx
		){

		float
		*d_UnpaddedKernel,
		*d_PaddedKernel;

		// Changed
		cufftHandle
		//	fftPlanFwd,
		//	fftPlanInv;
		fftPlan;

		bool bRetVal;

		#ifdef PRINT_KERNEL
			printf("Testing kernel spectrum computation\n");
			fprintf(stderr, "Testing kernel spectrum computation\n");
		#endif
		const int    kernelH = chem_cctx.KH;//kernel_cctx.DH;//7;
		const int    kernelW = chem_cctx.KW;//kernel_cctx.DW;//6;
		const int    kernelY = chem_cctx.KY;//kernel_cctx.DH / 2;//3;
		const int    kernelX = chem_cctx.KX;//kernel_cctx.DW / 2;//4;
		const int       fftH = chem_cctx.FFTH;
		const int       fftW = chem_cctx.FFTW;

	#ifdef PRINT_KERNEL
		printf("\tkernelH: %d\tkernelW: %d\n", kernelH, kernelW);
		printf("\tkernelX: %d\tkernelY: %d\n", kernelX, kernelY);
		printf("\tfftH: %d\tfftW: %d\n", fftH, fftW);
		fprintf(stderr,"...allocating memory\n");
	#endif
		checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_UnpaddedKernel, kernelH * kernelW * sizeof(float)));

		// Changed
	#ifdef PRINT_KERNEL
		//	printf("...creating R2C FFT plans for %i x %i\n", fftH, fftW);
	#endif
			//	checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
			//	checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));
	#ifdef PRINT_KERNEL
			printf("...creating C2C FFT plan for %i x %i\n", fftH, fftW / 2);
	#endif
		checkCudaErrors(cufftPlan2d(&fftPlan, fftH, fftW / 2, CUFFT_C2C));

	#ifdef PRINT_KERNEL
		fprintf(stderr,"...uploading to GPU and padding convolution kernel and input data\n");
	#endif
		checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fftH * fftW * sizeof(float)));
		padKernel(
				d_PaddedKernel,
				d_Kernel,
				fftH,
				fftW,
				kernelH,
				kernelW,
				kernelY,
				kernelX
		);

		checkCudaErrors(cudaDeviceSynchronize());

		// Changed
		// d_KernelSpectrum = FFT(d_PaddedKernel)
		//	checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel,       (cufftComplex *)d_KernelSpectrum));
		//CUFFT_INVERSE works just as well...
			const int FFT_DIR = CUFFT_FORWARD;
	#ifdef PRINT_KERNEL
		printf("...transforming convolution kernel\n");
	#endif
		checkCudaErrors(cufftExecC2C(fftPlan, (cufftComplex *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum, FFT_DIR));
		checkCudaErrors(cudaDeviceSynchronize());

		// Changed
		//	checkCudaErrors(cufftDestroy(fftPlanFwd));
		checkCudaErrors(cufftDestroy(fftPlan));
		checkCudaErrors(cudaFree(d_PaddedKernel));
		checkCudaErrors(cudaFree(d_UnpaddedKernel));

		return bRetVal;
	}


	bool computeChemDiffusionGPU(
			float		*h_ChemOut,
			float		*h_ChemIn,
	//		fComplex	*d_KernelSpectrum,
			fComplex	*d_KernelSpectrum0,
			c_ctx		 cctx,
			int		 	 iter)
	{
		float t = 0.0;

		float
		*d_Data,
		*d_PaddedData,
		*d_UnpaddedResult;
		//	*d_PaddedKernel;

		// Changed
		fComplex
	//	*d_DataSpectrum;
			*d_DataSpectrum0;

		// Changed
		cufftHandle
	//	fftPlanFwd,
		//	fftPlanInv;
		fftPlan;

		bool bRetVal;
		StopWatchInterface *hTimer = NULL;
		sdkCreateTimer(&hTimer);

	#ifdef PRINT_KERNEL
		printf("Testing GPU chemical diffusion computation\n");
		fprintf(stderr,"Testing GPU chemical diffusion computation\n");
	#endif
		const int    kernelH = cctx.KH;
		const int    kernelW = cctx.KW;
		const int    kernelY = cctx.KY;
		const int    kernelX = cctx.KX;
		const int      dataH = cctx.DH;
		const int      dataW = cctx.DW;
		const int       fftH = cctx.FFTH;
		const int       fftW = cctx.FFTW;

	#ifdef PRINT_KERNEL
		printf("\tkernelH: %d\tkernelW: %d\n", kernelH, kernelW);
		printf("\tkernelX: %d\tkernelY: %d\n", kernelX, kernelY);
		printf("\tdataH: %d\tdataW: %d\n", dataH, dataH);
		printf("\tfftH: %d\tfftW: %d\n", fftH, fftW);
		fprintf(stderr,"...allocating memory\n");
	#endif
		checkCudaErrors(cudaMalloc((void **)&d_Data,   dataH   * dataW   * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fftH * fftW * sizeof(float)));
		checkCudaErrors(cudaMalloc((void **)&d_UnpaddedResult, dataH * dataW * sizeof(float)));

		// Changed
		//	checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum,   fftH * (fftW / 2 + 1) * sizeof(fComplex)));
		checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum0,   fftH * (fftW / 2) * sizeof(fComplex)));

		// Changed
	#ifdef PRINT_KERNEL
		//	printf("...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
	#endif
		//	checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
		//	checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));
	#ifdef PRINT_KERNEL
		printf("...creating C2C FFT plan for %i x %i\n", fftH, fftW / 2);
	#endif
		checkCudaErrors(cufftPlan2d(&fftPlan, fftH, fftW / 2, CUFFT_C2C));

	#ifdef PRINT_KERNEL
		fprintf(stderr,"...uploading to GPU and padding input data\n");
	#endif
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		checkCudaErrors(cudaMemcpy(d_Data,   h_ChemIn,   dataH   * dataW *   sizeof(float), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemset(d_PaddedData,   0, fftH * fftW * sizeof(float)));
		sdkStopTimer(&hTimer);
		double dataTransferTime = sdkGetTimerValue(&hTimer);

	/*	padDataRightWall(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX);
	*/

	//	padDataMirror(
	//			d_PaddedData,
	//			d_Data,
	//			fftH,
	//			fftW,
	//			dataH,
	//			dataW,
	//			kernelH,
	//			kernelW,
	//			kernelY,
	//			kernelX
	//	);

	//	padDataClampToBorder(
	//			d_PaddedData,
	//			d_Data,
	//			fftH,
	//			fftW,
	//			dataH,
	//			dataW,
	//			kernelH,
	//			kernelW,
	//			kernelY,
	//			kernelX
	//	);

	#ifdef PRINT_KERNEL
		fprintf(stderr,"...performing convolution\n");
	#endif

		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);

		// Changed : Added
		//CUFFT_INVERSE works just as well...
		const int FFT_DIR = CUFFT_FORWARD;

		//	checkCudaErrors(cudaDeviceSynchronize());
		//	checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));
		checkCudaErrors(cudaDeviceSynchronize());

		// --------- Computing convolution ------------ begin

		// d_DataSpectrum = FFT(d_PaddedData)
		//	checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));
		checkCudaErrors(cufftExecC2C(fftPlan, (cufftComplex *)d_PaddedData, (cufftComplex *)d_DataSpectrum0, FFT_DIR));

		// d_DataSpectrum = d_DataSpectrum * d_KernelSpectrum
		//	modulateAndNormalize(d_DataSpectrum, d_KernelSpectrum, fftH, fftW, 1);
		#ifdef PRINT_KERNEL
			printf( "fftH: %d\tfftW: %d\n", fftH, fftW);
		#endif
		spProcess2D(d_DataSpectrum0, d_DataSpectrum0, d_KernelSpectrum0, fftH, fftW / 2, FFT_DIR);

		// d_PaddedData = IFFT(d_DataSpectrum)			<------- Output
		//	checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));
		checkCudaErrors(cufftExecC2C(fftPlan, (cufftComplex *)d_DataSpectrum0, (cufftComplex *)d_PaddedData, -FFT_DIR));
		// --------- Computing convolution ------------ end

		checkCudaErrors(cudaDeviceSynchronize());
		sdkStopTimer(&hTimer);
		double gpuTime = sdkGetTimerValue(&hTimer);
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);

	#ifdef PRINT_KERNEL
		fprintf(stderr,"...removing results padding\n");
	#endif
	/*	unpadResult(
				d_UnpaddedResult,
				d_PaddedData,
				dataH,
				dataW,
				fftH,
				fftW
		);
	*/
		sdkStopTimer(&hTimer);
		double unpadTime = sdkGetTimerValue(&hTimer);

	#ifdef PRINT_KERNEL
		fprintf(stderr,"...reading back GPU convolution results\n");
	#endif
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);
		// h_ChemOut = d_UnpaddedResult
		checkCudaErrors(cudaMemcpy(h_ChemOut, d_UnpaddedResult, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));
		sdkStopTimer(&hTimer);
		double readbackTime = sdkGetTimerValue(&hTimer);

	#ifdef PRINT_KERNEL
		printf("\t\tdata transfer:				%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
		printf("\t\tGPU chemical diffusion computation:	%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (gpuTime * 0.001), gpuTime);
		printf("\t\textract results:			%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (unpadTime * 0.001), unpadTime);
		printf("\t\tread back:				%f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (readbackTime * 0.001), readbackTime);
		fprintf(stderr,"...shutting down\n");
	#endif

		sdkDeleteTimer(&hTimer);

		// Changed
		//	checkCudaErrors(cufftDestroy(fftPlanInv));
		//	checkCudaErrors(cufftDestroy(fftPlanFwd));
		checkCudaErrors(cufftDestroy(fftPlan));
	#ifdef PRINT_KERNEL
		printf("...freeing device pointers\n");
	#endif
		// Changed
		//	checkCudaErrors(cudaFree(d_DataSpectrum));
		checkCudaErrors(cudaFree(d_DataSpectrum0));
		checkCudaErrors(cudaFree(d_PaddedData));
		checkCudaErrors(cudaFree(d_Data));
		checkCudaErrors(cudaFree(d_UnpaddedResult));

		//	free(h_ResultGPU);				// Use outside of function as input for next iteration
		//	free(h_Data);

	#ifdef PRINT_KERNEL
		printf("...returning to main()\n");
	#endif
		return bRetVal;
	}
#endif	// GPU_DIFFUSE (*)

using namespace std;

/* ------------------------------------------------------------------------------------ */
/*                            STATIC VARIABLES INITIALIZATIONS                          */
/* ------------------------------------------------------------------------------------ */
double WHWorld::clock = 0;
unsigned WHWorld::seed = 27000; //initial number of cells
bool WHWorld::highTNFdamage = false;
float WHWorld::patchpermm = 0;
#ifdef MODEL_SCAFFOLD
	int WHWorld::initialCaAlg = 0;
	float G;        // Elastic Modulus (kPa)
	float pXL;     // Crosslink Density (mmol/mL = M)
	float Alg_Mn ;  // molecular weight of alginate (kDa)
	float Q;        // Swelling Ratio
	float w = 0;    // Mass Loss (%)
	float poreWidth = 200.00;    // (um)
#endif

#ifdef MODEL_SCAFFOLD
	float WHWorld::Ca_Mw = 3400;      // Ca Molecular Weight (Mw ≈ 3,400 = g/mol)
	float WHWorld::Alg_Mn = 1500;     // Average molecular weight (Mw = 1 kDa = 1000 g/mol)
	//float WHWorld::Alg_Mn = 90;
	//float WHWorld::Alg_Mn = 200; //143;
#endif

float WHWorld::thresholdTNFdamage = 10.0; //ng
float WHWorld::cytokineDecay[6] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.5}; // 0.2, 0.2,
float WHWorld::halfLifes_static[6] = {33.6, 2.7, 46, 103, 24, 60}; // 13, 13,

#ifdef MODEL_SCAFFOLD
	float WHWorld::ElasticMod[7] = {125, 58, 971, 1.037, 756, 0.516, 0.165}; 
	float WHWorld::XLDensity[2] = {2.3, 10.1};
	float WHWorld::SwellRatio[5] = {72.478, 0.131, 22.034, 3.284, 35.752};	//float WHWorld::SwellRatio[5] = {0.4, 0.4, 3, 7.9, 1400};
	float WHWorld::MassLoss[4] = {0.234, 7.785, 0.15, 1.36};	//float WHWorld::MassLoss[4] = {17.6, 0.9, 60, 5.3};
	float WHWorld::PoreSize[2] = {1769.8, 258.5};	//float WHWorld::PoreSize[3] = {345.2, 309.9, 138.1};
#endif

WHWorld::WHWorld(double length, double width, double height, double plength) {
	// Generate random seeds:
	for(int i = 0; i < NUM_THREAD; i++) seeds[i] = 25234 + 17*i;

	// Allocate memory for local lists of cell pointers to add:
	for (int i = 0; i < MAX_NUM_THREADS; i++) localNewCells[i] = new vector<Cell*>;

	/* -------------------------------------------------------------------------- */
	/*                                 GRID SETUP                                 */
	/* -------------------------------------------------------------------------- */
	this->patchlength = plength;

    // Number of patches in x,y,z dimensions:
	int nx = width/patchlength;
	int ny = length/patchlength;
	int nz = height/patchlength;
	World::setupGrid(
        nx,             // number of grid points (patches) in x dimension
        ny,             // number of grid points (patches) in y dimension
        nz,             // number of grid points (patches) in z dimension
        0.0,            // min coordinates in x
        width,          // max corodinates in x
        0.0,            // min coordinates in y
        length,         // max coordinates in y
        0.0,            // min coordinates in z
        height          // max coordinates in z
    );

	// Read input parameters (chem baseline, wound dimensions, initial cells) from config file
	int temp = WHWorld::userInput();
	cout << "length, width, height: " << length << " mm, " << width << " mm, " << height << " mm" << endl;
	cout << "Number of patches: nx, ny, nz: " << nx << ", " << ny << ", " << nz << " " << endl;

	// Allocate and initialize Patches/ECM
	if (util::ABMerror(!(this->worldPatch = new Patch [(nx)*(ny)*(nz)]), "Patch mem alloc error!", __FILE__, __LINE__)) exit(1);
	if (util::ABMerror(!(this->worldECM = new ECM [(nx)*(ny)*(nz)]), "ECM mem alloc error!", __FILE__, __LINE__)) exit(1);
	cout << "worldPatch size: " << (nx)*(ny)*(nz) << " (Number of world patches) "<< endl;
	cout << "worldECM size: " << (nx)*(ny)*(nz) << " (Number of ECM patches) " << endl;

 	/* Try initializing Patches and ECMs with the threads that will access them later since the default allocation policy on Linux platforms is first-touch. 
	 * This is a best-effort implementation, since we cannot guarantee size of data accessed per thread to be an integer multiple of page size. */
    std::cout << std::fixed;
    std::cout << std::setprecision(3);
	//cout << "	allocating ECM Managers (also Patches) with best-effort first touch" << endl;
	
	for (int iz = 0; iz < nz; iz++) {
		#pragma omp parallel for
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
				this->worldPatch[in] = Patch(ix, iy, iz, in);
				this->worldECM[in] = ECM(ix, iy, iz, in);
			}
		}
	}

	// Define Class static variables and pointers
	WHWorld::patchpermm = nx/width;
	Agent::nx = this->nx;
	Agent::ny = this->ny;
	Agent::nz = this->nz;
	Agent::agentWorldPtr = this;
	Agent::agentPatchPtr = this->worldPatch;
	Agent::agentECMPtr = this->worldECM;
	ECM::ECMWorldPtr = this;
	ECM::ECMPatchPtr = this->worldPatch;

	/* ----------------------- INITIALIZATION SUBROUTINES ----------------------- */

    /* Define initial attributes of patches, damage, ECM, chem and cells based on user defined values (in config file) and traits of native tissue */
	this->initializePatches();
	this->initializeECM();
	this->initializeChem();
	this->initializeCells();
	#ifdef MODEL_SCAFFOLD
		this->initializeCaAlg();
		/* Create a temp Cell object to be able to call the Agent function cellCaAlgBehavior, as Agent is an abstract class */
		Cell tmpAgent;
		Cell* tmpThis = &tmpAgent;
		tmpThis->Agent::cellCaAlgBehavior();
		//Agent::cellCaAlgBehavior(); 
	#endif
	this->initializeDamage();

	/* Calling update functions to synchronize read and write portion of the attributes */
	//this->updateChemCPU();
	this->updateCellsInitial();  // Add cells to list before removal and updates
	this->updateECMManagers();
	this->updatePatches();
	//cout << "setupGrid complete" << endl;
    cout << "-------------------------------------------" << endl; 
}

WHWorld:: ~WHWorld(){
	free(this->D);
	free(this->HalfLifes);

	#ifdef GPU_DIFFUSE
		delete[] h_diffusion_results;
		free(this->chem_cctx);
	#endif
	for (int i = 0; i < MAX_NUM_THREADS; i++) delete localNewCells[i];
	cerr << " removing dead cells" << endl;
	
	int cellsSize = cells.size();
	
	#pragma omp parallel for
	for (int i = 0; i < cellsSize; i++) {
		#ifdef _OMP
			int tid = omp_get_thread_num();
		#else
			int tid = DEFAULT_TID;
		#endif

		Cell* cell = cells.getDataAt(i);
		if (!cell) continue;
		cells.deleteData(i, tid);
		delete cell;
 	}

	for (int ic = 0; ic < this->typesOfChem; ic++) if (chemAllocation[ic] != NULL) delete [] chemAllocation[ic];
    
    if (chemAllocation != NULL) delete [] chemAllocation;
    if (worldPatch != NULL)	delete [] worldPatch;
    if (worldECM != NULL) delete [] worldECM;
    cout << "WHWorld has been successfully destructed." << endl;
}

void destroyCell(Cell* &agent) {
	if (agent) {
		free(agent);
		agent = NULL;
	}
}

void WHWorld::assignPatches(int type, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax){
	// Assign patches within bounds of type 'type'
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
                switch (type){
                    case damage:
                        this->worldPatch[in].type[read_t] = damage;
						this->worldPatch[in].type[write_t] = damage;
                		this->worldPatch[in].color[read_t] = cdamage;
                		this->worldPatch[in].color[write_t] = cdamage;
                		this->worldPatch[in].dirty = true; 
                    	break; 
                    case CaAlg:
						this->worldPatch[in].type[read_t] = CaAlg;
						this->worldPatch[in].type[write_t] = CaAlg;
						this->worldPatch[in].color[read_t] = cCaAlg;
						this->worldPatch[in].color[write_t] = cCaAlg;
						this->worldPatch[in].dirty = true; 
						break; 
                }
            }
        }
    }
}

void WHWorld::initializePatches() {
  #ifdef MODEL_3D
        assignPatches(CaAlg, 0, nx, 0, ny, 0, nz); 
  #else
        assignPatches(CaAlg, 0, nx, 0, ny, 0, 0); 
  #endif

	// Assign values to initial:
	WHWorld::initialCaAlg = this->countPatchType(CaAlg);
	//cout << "Finished building Ca-Alg Hydrogel" << endl;
}

#ifdef MODEL_SCAFFOLD
	void WHWorld::initializeCaAlg(){
		cout << "Begin Calculating Ca-Alg Properties..." << endl; 
			
		/* ---------------------- Parameters of Ca-Alg Scaffold --------------------- */
		float Alg_ww = this->Alg_wv/(this->Alg_wv); 

		/* ----- Calculate initial bulk mechanical properties of Ca-Alg Scaffold ---- */
		/* p_XL: Crosslink Density (mmol/mL = M) 
		* 		 Linear dependence of Shear modulus on cross-link concentration for constant polymer concentration 
		*/

		this->pXL = 0.014;
		//this->pXL = 0.022;
		//this->pXL = 0.029;
		//this->pXL = 0.034;
		//this->pXL = 0.02;

		cout << "		Final Alginate concentration (%w/v): " << this->Alg_wv<< endl; 
		cout << "       Alginate Molecular Weight (kDa) = " << this->Alg_Mn << endl;
		cout << "       Calcium Crosslinking Density (mmol/mL = M) = " << this->pXL << endl;

		/* Calculate Initial Elastic Modulus E (kPa)
		*  E = a (( b*TotalProtein(w/v) + c)* Alg(w/w) + d*TP(w/v)) + e*(f*Alg(w/v) + g)*XL(w/w)
		*
		*       Follows rule of mixtures where stiffness of mixture is weight average of components.
		*       Linear dependence of modulus on cross-link concentration for constant polymer concentration 
		*/			
		#ifdef CALIBRATION
			this->E = -WHWorld::ElasticMod[0] + WHWorld::ElasticMod[1]*(Alg_wv) - WHWorld::ElasticMod[2]*(pXL) + WHWorld::ElasticMod[3]*(Alg_Mn) + WHWorld::ElasticMod[4]*(Alg_wv)*(pXL) - WHWorld::ElasticMod[5]*(Alg_wv)*(Alg_Mn) - WHWorld::ElasticMod[6]*(pXL)*(Alg_Mn);
		#else 
			this->E = -125 + 58*(Alg_wv) - 971*(pXL) + 1.037*(Alg_Mn) + 756*(Alg_wv*pXL) - 0.516*(Alg_wv*Alg_Mn) - 0.165*(pXL*Alg_Mn);
		#endif
		cout << "       Elastic Modulus (kPa) = " << this->E << endl; 
		
		/* Pore Size (um): poreWidth = -a * Alg_ww^2 + b * Alg_ww + c */
		#ifdef CALIBRATION
			this->poreWidth = -WHWorld::PoreSize[0]*(pXL) + WHWorld::PoreSize[1]; 
		#else
			this->poreWidth = -1769.84*(pXL) + 258.5;  //-0.3113*pow(Alg_ww,2) + 1.5*Alg_ww + 50;   //this->poreWidth = (-0.01)*345.2*pow(Alg_ww,2) + 309.9*Alg_ww + 138.1; 
		#endif	
		cout << "     this->poreWidth = -" <<WHWorld::PoreSize[0]<<"*"<<(pXL)<<" + "<< WHWorld::PoreSize[1] << endl; 
		cout << "       Pore Width (um): " << this->poreWidth << endl; 
		
		/* Swelling Ratio:
		*       Swelling Ratio increase with Alg content and with time
		*       Important in retaining water, facilitating diffusion.
		*/
		#ifdef CALIBRATION
			this->Q = WHWorld::SwellRatio[0] - WHWorld::SwellRatio[1]*(this->reportDay()) - WHWorld::SwellRatio[2]*(Alg_wv)	- WHWorld::SwellRatio[3]*(this->reportDay())*(pXL) + WHWorld::SwellRatio[4]*(Alg_wv)*(pXL);
		#else
			this->Q = 72.478 - 0.131*(this->reportDay()) - 22.034*(Alg_wv) - 3.284*(this->reportDay())*(pXL) + 35.752*(Alg_wv)*(pXL);
		#endif 
		//cout << "    this->Q ="<< WHWorld::SwellRatio[0]<< "-"<< WHWorld::SwellRatio[1]<<"*"<<(this->reportDay())<<" - "<< WHWorld::SwellRatio[2]<<"*"<<(Alg_wv)<<" -"<< WHWorld::SwellRatio[3]<<"*"<<(this->reportDay())<<"*"<<(pXL) <<" + "<< WHWorld::SwellRatio[4]<<"*"<<(Alg_wv)<<"*"<<(pXL) << endl;
		cout << "       Swelling Ratio: " << this->Q << endl;

		/* Mass Loss (%) 
		*       Instable hydrogel degrades, replaced with cell-synthesized ECM proteins 
		*       Mass loss fraction (%) increases with time and Alg content
		*/
		#ifdef CALIBRATION
			this->w = 0;// this->w = WHWorld::MassLoss[0] + WHWorld::MassLoss[1]*(pXL) + WHWorld::MassLoss[2]*(reportDay()) - WHWorld::MassLoss[3]*(pXL)*(reportDay());
		#else
			this->w = 0.234 + 7.785*(pXL) + 0.15*(reportDay()) - 1.36*(pXL)*(reportDay());
		#endif
		if (w < 0 ) w = 0;  // no negative mass loss
		//cout << " 	this->w ="<< WHWorld::MassLoss[0]<< " +"<< WHWorld::MassLoss[1]<<"*"<<(pXL)<<" +"<< WHWorld::MassLoss[2]<<"*"<<(reportDay())<<" - "<< WHWorld::MassLoss[3]<<"*"<<(pXL)<<"*"<<(reportDay()) << endl; 
		cout << " Mass Loss (%): " << this->w << endl; 

	cout << "Finished calculating initial Ca-Alg properties" << endl;
	}
#endif //MODEL_SCAFFOLD

void WHWorld::initializeChem(){
	#ifdef GPU_DIFFUSE
		this->initializeChemGPU();
	#else
		this->initializeChemCPU();
	#endif
}

void WHWorld::initializeChemCPU() {
	/* Allocate two-dimensional matrix (chemAllocation[chemtype][patch index]) to store quantities of each chemical at each patch */
	this->typesOfChem = (this->baselineChem.size())*2 + 3;
	this->chemAllocation = new float*[this->typesOfChem];
	//cout << "number of chem allocated is "<< this->typesOfChem << endl;
		
	for (int ic = 0; ic < this->typesOfChem; ic++) {
		if (util::ABMerror(!(this->chemAllocation[ic] = new float[nx*ny*nz] ), "InitializeChem mem alloc error!", __FILE__, __LINE__)) exit(1);
	}

	// Link World attribute chemAllocation with WHChemical class attribute WHWorldChem:
	this->WHWorldChem.pTNF = this->chemAllocation[pTNF];
	this->WHWorldChem.pTGF = this->chemAllocation[pTGF];
	this->WHWorldChem.pIL1beta = this->chemAllocation[pIL1beta];
	this->WHWorldChem.dTNF = this->chemAllocation[dTNF];
	this->WHWorldChem.dTGF = this->chemAllocation[dTGF];
	this->WHWorldChem.dIL1beta = this->chemAllocation[dIL1beta];
	this->WHWorldChem.pcellgrad = this->chemAllocation[pcellgrad];

	// Initialize chemical concentrations:
	this->WHWorldChem.totalTNF = 0;
	this->WHWorldChem.totalTGF = 0;
	this->WHWorldChem.totalIL1beta = 0;
		
	int countCaAlg = this->countPatchType(CaAlg);
	if (this->baselineChem.size() == 4) {
		for (int iz = 0; iz < this->nz; iz++) {
			/* Try initializing chemicals with the threads that will access them later since the default allocation policy on Linux platforms is first-touch.
		 	This is a best-effort implementation, since we cannot guarantee size of data accessed per thread to be an integer multiple of page size. */
			#pragma omp parallel for
			for (int iy = 0; iy < this->ny; iy++) {
				for (int ix = 0; ix < this->nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;
					this->WHWorldChem.dTNF[in] = 0;
					this->WHWorldChem.dTGF[in] = 0;
					this->WHWorldChem.dIL1beta[in] = 0;

					// Baseline chemical concentrations are initialized in tissue 
					if (this->worldPatch[in].type[read_t] == CaAlg) {
						this->WHWorldChem.pTNF[in] = this->baselineChem[TNF]/countCaAlg;
						this->WHWorldChem.pTGF[in] = this->baselineChem[TGF]/countCaAlg;
						this->WHWorldChem.pIL1beta[in] = this->baselineChem[IL1beta]/countCaAlg;
					} else {
						this->WHWorldChem.pTNF[in] = 0;
						this->WHWorldChem.pTGF[in] = 0;
						this->WHWorldChem.pIL1beta[in] = 0;
					}
						
					// Initialize chemical gradient levels that agents are attracted by 
					float patchIL1 = this->WHWorldChem.pIL1beta[in];
					float patchTNF = this->WHWorldChem.pTNF[in];
					float patchTGF = this->WHWorldChem.pTGF[in];
					float patchcollagen = this->worldECM[in].fcollagen[read_t];
					//float grad = patchTNF + patchTGF + patchcollagen + patchIL1;
					this->WHWorldChem.pcellgrad[in] = patchTGF;  
					#pragma omp critical
					{
						//Initialize total chemical concentration:
						this->WHWorldChem.totalTNF += this->WHWorldChem.pTNF[in];
						this->WHWorldChem.totalTGF += this->WHWorldChem.pTGF[in];
						this->WHWorldChem.totalIL1beta += this->WHWorldChem.pIL1beta[in];
					}
				}
			}
		}
		cout << "		Initial cytokine concentrations: totalTNF = " << this->WHWorldChem.totalTNF << ", totalTGF = " << this->WHWorldChem.totalTGF << ", totalIL1beta = " << this->WHWorldChem.totalIL1beta << endl;
	} else if (util::ABMerror(1, "Error initializing chemicals!!", __FILE__, __LINE__)) exit(1);
	//cout << "Finished initializing chem" << endl;
}

void printWindow(float* a, int h, int w, int r){
	#ifdef PRINT_KERNEL
		int mx = w/2;
		int my = h/2;
		int i_start	= my - r;
		int i_end	= my + r + 1;
		int j_start	= mx - r;
		int j_end	= mx + r + 1;
		if (i_start < 0 || i_end >= h || j_start < 0 || j_end >= w) {
			cout << "printWindowError: " << h << "x" << w << endl;
			cout << "	i_start:	" << i_start << endl;
			cout << "	i_end:		" << i_end << endl;
			cout << "	j_start:	" << j_start << endl;
			cout << "	j_end:		" << j_end << endl;
			exit(-1);
		}
		cout << "Matrix window: ---" << endl;
		for (int i = i_start; i < i_end; i++) {
			for (int j = j_start; j < j_end; j++) cout << a[i*w + j] << " ";
			cout << endl;
		}
	#endif
}

#ifdef GPU_DIFFUSE
	void WHWorld::initializeChemGPU() {
		/********************************************
		* Kernel Computation	                    *
		********************************************/
		int numChem = baselineChem.size();
		// Prepare dimensions
		int nSec = 1800;//90;//180;//600;
		// TODO: calculate dt w.r.t. max D (for now assume max D = 20.0)
		float dt = 2.5;
		float dx = 15.0;	// um
		//float D = 20.0;
		float dx2 = dx * dx;
		//float lambda = (D*dt)/dx2;
		float *lambda = (float *) malloc(numChem * sizeof(float));
		float *gamma = (float *) malloc(numChem * sizeof(float));

		#ifdef PRINT_KERNEL
			cerr << "Allocated " << numChem << " elements for lambda (" << lambda << ")  and gamma (" << gamma << ")" << endl;
			cerr << "D: " << D << "		halflife: " << HalfLifes << endl;
		#endif

		for (int ic = 0; ic < numChem; ic++) {
			lambda[ic] = (D[ic]*dt)/dx2;
			// Overwrite halflifes from config file. Use instead values from sensitivity analysis input file (Sample.txt)
			HalfLifes[ic] = halfLifes_static[ic] * 60;		// minutes * 60
			gamma[ic] = 1 - pow(2, -(1/HalfLifes[ic]));
		}

		#ifdef PRINT_KERNEL
			cerr << "Done--" << endl;
		#endif

		int kernelRadius	= nSec/dt;
		const int outKernelH		= (kernelRadius-1) * 2 + 1;
		const int outKernelW		= outKernelH;
		const int    kfftH = snapTransformSize(outKernelH + KCoeffsH - 1);
		const int    kfftW = snapTransformSize(outKernelW + KCoeffsW - 1);

		// Changed 2 : Added
		//	const int windowH = pow(2, ceil(log(outKernelH/2)/log(2)));
		//	const int windowW = pow(2, ceil(log(outKernelW/2)/log(2)));
		const int windowH = (kernelRadius/2) * 2 + 1; //outKernelH/2;
		const int windowW = (kernelRadius/2) * 2 + 1; //outKernelW/2;

		// Kernel convolution context
		c_ctx kernel_cctx;

		// Set dimensions
		kernel_cctx.DH		= outKernelH;
		kernel_cctx.DW		= outKernelW;
		kernel_cctx.KH		= KCoeffsH;
		kernel_cctx.KW		= KCoeffsW;
		kernel_cctx.KX		= 1;//0;
		kernel_cctx.KY		= 1;//0;
		kernel_cctx.FFTH	= kfftH;
		kernel_cctx.FFTW	= kfftW;
		// Changed 2
		kernel_cctx.windowH = windowH;
		kernel_cctx.windowW = windowW;
		//float gamma		= 0.00019252;//0.00;		// half life = 60 minutes

		float
		*h_Kernel,						// In place of h_ResultGPU
		**d_Kernel,						// In place of d_UnpaddedResult -> final kernel
		**h_dKernel;

		// Changed 2
		float
		*h_Window,
		**d_Window,
		**h_dWindow;

		//	// Changed
		//	fComplex
		////	*d_KernelSpectrum;
		// 	*d_KernelSpectrum0;

		// Changed 2
		h_Window = (float *)malloc(windowH * windowW    * sizeof(float));
		// Allocate memory for kernel computation
		h_Kernel	= (float *)malloc(outKernelH    * outKernelW * sizeof(float));
		// Allocate array of 'numChem' pointers for device kernels
		checkCudaErrors(cudaMalloc((void **)&d_Kernel, numChem * sizeof(float*)));
		h_dKernel = (float **) malloc (numChem * sizeof(float*));
		checkCudaErrors(cudaMemcpy(h_dKernel,d_Kernel,numChem * sizeof(float*),cudaMemcpyDeviceToHost));
		// Changed 2 : Added
		// Allocate array of 'numChem' pointers for device kernel center windows
		checkCudaErrors(cudaMalloc((void **)&d_Window, numChem * sizeof(float*)));
		h_dWindow = (float **) malloc (numChem * sizeof(float*));
		checkCudaErrors(cudaMemcpy(h_dWindow,d_Window,numChem * sizeof(float*),cudaMemcpyDeviceToHost));

		// Allocate 'numChem' arrays for device kernels
		for (int ic = 0; ic < numChem; ic++) {
			#ifdef PRINT_KERNEL
				cout << "Allocating kernel arrays" << endl;
			#endif
			checkCudaErrors(cudaMalloc((void **)&h_dKernel[ic], outKernelH * outKernelW * sizeof(float)));
			checkCudaErrors(cudaMalloc((void **)&h_dWindow[ic], windowH * windowW * sizeof(float)));
			#ifdef PRINT_KERNEL
				cout << "Computing kernel arrays " << ic << endl;
			#endif
			// KERNEL COMPUTATION
			if (!computeKernel(h_Kernel,h_dKernel[ic],h_Window,h_dWindow[ic],kernelRadius,lambda[ic],gamma[ic],dt,kernel_cctx)){} // TODO: Error Handling
		}
		#ifdef PRINT_KERNEL
		cout << "Done---" << endl;
		#endif

		printWindow(h_Kernel, outKernelH, outKernelW, 4);
		printWindow(h_Window, windowH, windowW, 4);

		// Free pointers
		if (h_Kernel) free(h_Kernel);	h_Kernel = NULL;
		if (h_Window) free(h_Window);	h_Window = NULL;

		/********************************************
		* Kernel Spectrum Computation	            *
		********************************************/
		const int    chemH = this->ny;// - 1;
		const int    chemW = this->nx;// - 1;

		// Changed 2
		//	const int    cfftH = snapTransformSize(chemH + kernel_cctx.FFTH - 1);
		//	const int    cfftW = snapTransformSize(chemW + kernel_cctx.FFTW - 1);
		const int    cfftH = snapTransformSize(chemH + kernel_cctx.windowH - 1);
		const int    cfftW = snapTransformSize(chemW + kernel_cctx.windowW - 1);
		//printf("Chem world: %d x %d\n", chemH, chemW);
		//fComplex** h_dKernel_spectrum;

		// Allocate memory on GPU for diffusion kernel spectrums for convolution
		checkCudaErrors(cudaMalloc((void **)&d_kernel_spectrum,numChem * sizeof(fComplex *)));
		h_dKernel_spectrum = (fComplex **) malloc(numChem * sizeof(fComplex*));
		checkCudaErrors(cudaMemcpy(this->h_dKernel_spectrum,d_kernel_spectrum,numChem * sizeof(fComplex*),cudaMemcpyDeviceToHost));

		for (int ic = 0; ic < numChem; ic++) {
			#ifdef PRINT_KERNEL
				cout << "Allocation: " << ic << endl;
			#endif
			checkCudaErrors(cudaMalloc(
				(void **)&(this->h_dKernel_spectrum[ic]),
				// Changed
				//cfftH * (cfftW / 2 + 1) * sizeof(fComplex))
				cfftH * (cfftW / 2) * sizeof(fComplex))
			);
		}
		#ifdef PRINT_KERNEL
			cout << "Done" << endl;
		#endif

		// Allocate and set dimensions in chemical convolution context
		this->chem_cctx         = (c_ctx*) malloc(sizeof(c_ctx));
		this->chem_cctx->DH	= chemH;
		this->chem_cctx->DW	= chemW;

		// Changed 2
		this->chem_cctx->KH	= kernel_cctx.windowH;
		this->chem_cctx->KW	= kernel_cctx.windowW;
		this->chem_cctx->KX	= kernel_cctx.windowW / 2;
		this->chem_cctx->KY	= kernel_cctx.windowH / 2;
		//	this->chem_cctx->KH	= kernel_cctx.DH;//kernel_cctx.FFTH;
		//	this->chem_cctx->KW	= kernel_cctx.DW;//kernel_cctx.FFTW;
		//	this->chem_cctx->KX	= kernelRadius-1;//128;//1;
		//	this->chem_cctx->KY	= kernelRadius-1;//128;//1;

		this->chem_cctx->FFTH	= cfftH;
		this->chem_cctx->FFTW	= cfftW;
		this->chem_cctx->windowH = -1;
		this->chem_cctx->windowW = -1;

		StopWatchInterface *hTimer = NULL;
		sdkCreateTimer(&hTimer);
		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);

		// COMPUTE KERNEL SPECTRUM
		for (int ic = 0; ic < numChem; ic++) {
			computeKernelSpectrum(
				this->h_dKernel_spectrum[ic],
				// Changed 2
				//h_dKernel[ic],
				h_dWindow[ic],
				kernel_cctx,
				*(this->chem_cctx)
			);
			// d_Kernel not used anymore, only need kernel spectrum
			checkCudaErrors(cudaFree(h_dKernel[ic]));
		}

		sdkStopTimer(&hTimer);
		double kernelSpectrumComputationTime = sdkGetTimerValue(&hTimer);
		//printf("\tTotal kernel computation: %f MPix/s (%f ms)\n",(double)chemH * (double)chemW * 1e-6 / (kernelSpectrumComputationTime * 0.001),kernelSpectrumComputationTime);

		/********************************************
		* GPU Diffusion Preparation	            *
		********************************************/
		/* Allocate two-dimensional matrix (chemAllocation[chemtype][patch index]) to store quantities of each chemical at each patch */
		int base_chem_types       = this->baselineChem.size();
		this->typesOfChem         = base_chem_types*2 + 3;
		this->chemAllocation      = new float*[this->typesOfChem];

		// Allocate buffer to store diffusion results from GPU
		this->h_diffusion_results = new float*[base_chem_types];

		//cout << "number of chem allocated is "<< this->typesOfChem << endl;
		//cout << "number of chem allocated for GPU diffusion is " << base_chem_types << endl;
		for (int ic = 0; ic < this->typesOfChem; ic++) {
			if (util::ABMerror(!(this->chemAllocation[ic]  = new float[nx*ny*nz] ),"InitializeChem mem alloc error!",__FILE__,__LINE__))exit(1);
			if (ic < base_chem_types)
				if (util::ABMerror(!(this->h_diffusion_results[ic] = new float[nx*ny*nz] ),"InitializeChem mem alloc error!",__FILE__,__LINE__))exit(1);
		}

		/* Link World attribute chemAllocation with WHWorldChem (WHWorldChem is linked to WHChemical) */
		this->WHWorldChem.pTNF = this->chemAllocation[pTNF];
		this->WHWorldChem.pTGF = this->chemAllocation[pTGF];
		this->WHWorldChem.pIL1beta = this->chemAllocation[pIL1beta];
		
		this->WHWorldChem.dTNF = this->chemAllocation[dTNF];
		this->WHWorldChem.dTGF = this->chemAllocation[dTGF];
		this->WHWorldChem.dIL1beta = this->chemAllocation[dIL1beta];
		this->WHWorldChem.pcellgrad = this->chemAllocation[pcellgrad];

		this->WHWorldChem.tTNF     = this->h_diffusion_results[pTNF];
		this->WHWorldChem.tTGF     = this->h_diffusion_results[pTGF];
		this->WHWorldChem.tIL1beta = this->h_diffusion_results[pIL1beta];

		// Initialize chemical concentrations:
		this->WHWorldChem.totalTNF = 0;
		this->WHWorldChem.totalTGF = 0;
		this->WHWorldChem.totalIL1beta = 0;
		
		int countCaAlg = WHWorld::initialCaAlg;

		if (this->baselineChem.size() == 3) {
			for (int iz = 0; iz < this->nz; iz++) {
			/* Try initializing chemicals with the threads that will access them later since the default allocation policy on Linux platforms is first-touch. 
			This is a best-effort implementation, since we cannot guarantee size of data accessed per thread to be an integer multiple  of page size. */
			#pragma omp parallel for
				for (int iy = 0; iy < this->ny; iy++) {
					for (int ix = 0; ix < this->nx; ix++) {
						int in = ix + iy*nx + iz*nx*ny;
						this->WHWorldChem.dTNF[in] = 0;
						this->WHWorldChem.dTGF[in] = 0;
						this->WHWorldChem.dIL1beta[in] = 0;

						// Baseline chemical concentrations are initialized in tissue:
						if (this->worldPatch[in].type[read_t] == CaAlg) {
							this->WHWorldChem.pTNF[in] = this->baselineChem[TNF]/countCaAlg;
							this->WHWorldChem.pTGF[in] = this->baselineChem[TGF]/countCaAlg;
							this->WHWorldChem.pIL1beta[in] = this->baselineChem[IL1beta]/countCaAlg;
						} else {
							this->WHWorldChem.pTNF[in] = 0;
							this->WHWorldChem.pTGF[in] = 0;
							this->WHWorldChem.pIL1beta[in] = 0;
						}

						// Initialize chemical gradient levels that agents are attracted by:
						float patchIL1 = this->WHWorldChem.pIL1beta[in];
						float patchTNF = this->WHWorldChem.pTNF[in];
						float patchTGF = this->WHWorldChem.pTGF[in];
						float patchcollagen = this->worldECM[in].fcollagen[read_t];
						//float grad = patchIL1 + patchTNF + patchTGF + patchFGF + patchcollagen;
						this->WHWorldChem.pcellgrad[in] = patchTGF;

						#pragma omp critical
						{
							//Initialize total chemical concentration:
							this->WHWorldChem.totalTNF += this->WHWorldChem.pTNF[in];
							this->WHWorldChem.totalTGF += this->WHWorldChem.pTGF[in];
							this->WHWorldChem.totalIL1beta += this->WHWorldChem.pIL1beta[in];
						}
					}
				}
			}
			//cout << "		Results from inside initialization:      totalTNF = " << this->WHWorldChem.totalTNF << ", totalTGF = " << this->WHWorldChem.totalTGF << ", totalIL1beta = " << this->WHWorldChem.totalIL1beta << endl;

		} else if (util::ABMerror(1,"Error initializing chemicals!!",__FILE__,__LINE__))exit(1);
		//cout << "Finished initializing chem" << endl;

		checkCudaErrors(cudaFree(d_Kernel));
		checkCudaErrors(cudaFree(d_Window));
		checkCudaErrors(cudaFree(d_kernel_spectrum));
		for (int ic = 0; ic < numChem; ic++) checkCudaErrors(cudaFree(h_dWindow[ic]));

		sdkDeleteTimer(&hTimer);
		free(lambda);
		free(gamma);
		free(h_Window);
		free(h_Kernel);
		free(h_dKernel);
		free(h_dWindow);

	}
#endif // GPU_DIFFUSE

void WHWorld::initializeCells() {
    //cout << "Begin Initializing Cells..." << endl;

	// Instantiate Cell list:
	cells = ArrayChain<Cell*>(DEFAULT_DATA_SMALL, 4, NULL, NULL);  // WHWorld::destroyChond);
    //cout << "Initialize Cells..." << endl; 

    // If initial cell count not input, seed scaffold with stem cells at density 10^6 cell/mL
    double cellDensity;     // cells/mm^3
	double scaffoldVolume = (nx*ny*nz)*pow(this->patchlength, 3);   //mm^3 
	double hydrogelVolume = scaffoldVolume*pow(10,-3); //mL

    if (this->initialCells[0] == 0){
        double cellpermL = 1*pow(10,6); 	 // (Xuan Li) 1 mill/ml
        cellDensity = cellpermL*pow(10,-3);  // cells/mm^3
    } else{
        double cellpermL = this->initialCells[0]/scaffoldVolume;
        cellDensity = cellpermL*pow(10,-3); 
    }
    
	this->initialCells[0] = cellDensity*scaffoldVolume; 
    int initialScaffoldCells = this->initialCells[0]; 
	
	cout << "		Scaffold Volume: " << scaffoldVolume << " mm^3 (" << hydrogelVolume << " mL)" << endl; 	
	cout << "		Cell Density: " << cellDensity << " cells/mm^3  (" << 1.0*pow(10,6) << " cells/mL)" << endl;
    cout << " 		Seeding " << initialScaffoldCells << " cells in scaffold " << endl; 

    // Sprout cell seeded hydrogel with mesenchymal stem cells at density 10^6 cells/mL
    sproutAgent(                    
        initialScaffoldCells,   // Number of cells to sprout
        CaAlg,                  // Type of patch to sprout on
        stem,            // Type of agent to sprout
        // Physical Boundary:
        0,                      //  -- left
        nx,                     //  -- right
        0,                      //  -- top
        ny,                     //  -- bottom
        0,                      //  -- front
        nz                      //  -- rear
	);
	//cout << "Finished initializing cells" << endl;
}

void WHWorld::initializeECM() {

	/* -------------------------------------------------------------------------- */
	/*                                  COLLAGEN                                  */
	/* -------------------------------------------------------------------------- */
   // cout << "Begin Initializing Collagen..." << endl;
    sproutAgentInArea(  nx*ny*nz,                     // Number of agents to sprout
            			CaAlg,                        // patch type
                        new_coll,                     // new collagen agent type
                        0,                            // lowest x-coordinate agent could be sprouted at
                        nx,                           // highest x-coordinate agent
                        0,                            // lowest y-coordinate
                        ny,                           // highest y-coordinate
                        0,                            // lowest z-coordinate
                        nz                            // highest z-coordinate
    );    // collagen in CaAlg Scaffold


	/* -------------------------------------------------------------------------- */
	/*                                  AGGRECAN                                  */
	/* -------------------------------------------------------------------------- */
    //cout << "Begin Initializing Aggrecan..." << endl;
    sproutAgentInArea(nx*ny*nz,                     // Number of agents to sprout
                      CaAlg,                        // patch type
                      new_agg,                      // new aggrecan agent type
                      0,                            // lowest x-coordinate agent could be sprouted at
                      nx,                           // highest x-coordinate agent
                      0,                            // lowest y-coordinate
                      ny,                           // highest y-coordinate
                      0,                            // lowest z-coordinate
                      nz                            // highest z-coordinate
    );    // aggrecan in CaAlg Scaffold
}

void WHWorld::initializeDamage() {
    // Sprout Damage throughout Scaffold and fragment ECM proteins 
    for (int iz = 0; iz < nz; iz++){
        for (int iy = 0; iy < ny; iy++){
            for (int ix = 0; ix < nx; ix++){
                int in = ix + iy*nx + iz*nx*ny; 
                worldPatch[in].inDamzone = true;
		        worldPatch[in].health[write_t] = 0;
		        worldPatch[in].damage[write_t] = worldPatch[in].damage[read_t]++;
                worldPatch[in].dirty = true;
				
				if (worldECM[in].empty[write_t] == true) continue;
				this->worldECM[in].fragmentNCollagen();
				this->worldECM[in].fragmentNAggrecan();
                worldPatch[in].updatePatch();
            }
        }
     }
	//cout <<  this->countPatchType(damage) << " damage created" << endl;
}

/*
 * Each call to WHWorld::go() performs the following major steps:
 * 	(0) Cell seedings
 * 	(1) Chemical diffusion
 * 	(2) Cell function
 * 	(3) ECM function
 * 	(4) Attributes synchronization
 * 			a) Update chemicals
 * 			b) Update cells
 * 			c) Update ECM managers
 * 			d) Update patches
*/
int WHWorld::go() {
    cout << "-------------------------------------------" << endl; 

	Patch* tempPatchPtr;
	Agent* tempAgentPtr;
	double hours = this->reportHour();
	double days = this->reportDay();

	// Profiling options defined in common.h
	#ifdef PROFILE_MAJOR_STEPS
		struct timeval start, end;
		long elapsed_time;  // in milliseconds
	#endif

	// Increment Clock in ticks (1 tick = 30 min)
	WHWorld::clock++;
	cout << "tick: " << clock << " , hour: " << hours << " , day: " << days << endl;

	#ifdef PROFILE_MAJOR_STEPS
		#if defined(GPU_DIFFUSE) && defined(_OMP)

		/* TIME_STAGE is a macro for timing a command/function and printing the timing info (See Utilities/time.h) */
		#pragma omp parallel num_threads(2)
		{
			int tid = omp_get_thread_num();
			if (tid == 1) {
				cout << "go() task 1" << endl;
				/* --------------------------- CHEMICAL DIFFUSION --------------------------- */
				TIME_STAGE(this->diffuseCytokines(), "Chemical diffusion", "1");
			
			} else if (tid == 0) {
				cout << "go() task 2" << endl;
				/* ------------------------------ CELL SEEDING ------------------------------ */
				//TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");

				/* ------------------------------ CELL FUNCTION ----------------------------- */
				TIME_STAGE(this->executeCells(), "Cells function", "2");

				/* ------------------------------ ECM FUNCTION ------------------------------ */
				TIME_STAGE(this->executeECMs(), "ECM function", "3(a)");
				TIME_STAGE(this->requestECMfragments(), "ECM fragmentation", "3(b)"); // Request fragment<ECM>

				/* ----------------------- ATTRIBUTES SYNCHRONIZATION ----------------------- */
				cerr << " begin update..." << endl;
				TIME_STAGE(this->updateCells(), "Update cells", "4(b)");
				TIME_STAGE(this->updateECMManagers(), "Update ECM Managers", "4(c)");
				TIME_STAGE(this->updatePatches(), "Update Patches", "4(d)");
			}
		}

		TIME_STAGE(this->updateChem(), "Update chem", "4(a)");

	#else	// GPU_DIFFUSE && _OMP
		/* TIME_STAGE is a macro for timing a command/function and printing the timing info (See Utilities/time.h) */

		/* --------------------------- CHEMICAL DIFFUSION --------------------------- */
		TIME_STAGE(this->diffuseCytokines(), "Chemical diffusion", "1");
		
		/* ------------------------------ CELL SEEDING ------------------------------ */
		//TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");
		
		/* ------------------------------ CELL FUNCTION ----------------------------- */
		TIME_STAGE(this->executeCells(), "Cells function", "2");

		/* ------------------------------ ECM FUNCTION ------------------------------ */
		TIME_STAGE(this->executeECMs(), "ECM function", "3(a)");
		TIME_STAGE(this->requestECMfragments(), "ECM fragmentation", "3(b)"); // Request fragment<ECM>

		/* ----------------------- ATTRIBUTES SYNCHRONIZATION ----------------------- */
		cerr << " begin update... " << endl;
		TIME_STAGE(this->updateCells(), "Update cells", "4(b)");
		TIME_STAGE(this->updateECMManagers(), "Update ECM Managers", "4(c)");
		TIME_STAGE(this->updatePatches(), "Update Patches", "4(d)");
		TIME_STAGE(this->updateChem(), "Update chem", "4(a)");
	#endif	// GPU_DIFFUSE && _OMP
	#else	// PROFILE_MAJOR_STEPS

		// For testing purposes:
		this->updateTotalChem();
		cout << "  TNF: " << this->WHWorldChem.totalTNF << ", TGF: " << this->WHWorldChem.totalTGF << ", IL1beta: " << this->WHWorldChem.totalIL1beta << endl;
		
		/* --------------------------- CHEMICAL DIFFUSION --------------------------- */
		this->diffuseCytokines();

		/* ------------------------------ CELL SEEDING ------------------------------ */
		//this->seedCells(hours);

		/* ------------------------------ CELL FUNCTION ----------------------------- */
		#ifdef MODEL_SCAFFOLD
			/* Create a temp Cell object to be able to call the Agent function cellCaAlgBehavior, as Agent is an abstract class */
			Cell tmpAgent;
			Cell* tmpThis = &tmpAgent;
			tmpThis->Agent::cellCaAlgBehavior();
			//Agent::cellCaAlgBehavior(); 
		#endif
		this->executeCells();

		/* ------------------------------ ECM FUNCTION ------------------------------ */
		this->executeECMs();
		this->requestECMfragments();

		#ifdef MODEL_SCAFFOLD
			/* ------------------------- UPDATE CaAlg Properties ------------------------ */
			this->updateSwellingRatio();
			this->updateMassLoss();
		#endif

		/* ----------------------- ATTRIBUTES SYNCHRONIZATION ----------------------- */
		cerr << " begin update... " << endl;
		this->updateCells();
		this->updateECMManagers();
		this->updatePatches();
		this->updateChem();

	#endif
	return 0;
}

/* -------------------------------------------------------------------------- */
/*                      MAJOR SECTION SUBROUTINES - begin                     */
/* -------------------------------------------------------------------------- */
void WHWorld::diffuseCytokines() {
	#ifdef PDE_DIFFUSE
	#ifdef GPU_DIFFUSE
		this->diffuseChemGPU();
	#else
		float D, deltat; 
			for(int ichem=0; ichem < this->baselineChem.size(); ichem++){
				/* Diffusion coefficient (mm^2/min) according to approximate Molecular weight. Examples of diffusion coefficients can be found at http://www.math.ubc.ca/~ais/website/status/diffuse.html
				* In Ca-Alg Scaffold, diffusion effected by hydrogel swelling (Q %) */
				switch(ichem){  // cytokine specific diffusion coefficient
					case TNF:                                                                                                               
						D = 0.0018*0.1*this->Q;                 
					case TGF:
						D = 0.00156*0.1*this->Q;                
					case IL1beta:
						D = 0.0018*0.1*this->Q;
					/* case TNF:                                                                                                               
						D = 0.0018;
					case TGF:
						D = 0.00156;                
					case IL1beta:
						D = 0.0018;
					*/
				}

				#ifdef MODEL_3D
					// 3D diffusion stability condition, dt < dx^2/6*D min given patchlength (dx = dy = dz), D from above
					deltat = pow(this->patchlength,2)/(6*D)-0.001; 
				#else
					// 2D diffusion stability condition, dt < dx^2/4*D min = 0.03125 assuming dx = dy = 0.015 mm, D from above
					deltat = pow(this->patchlength,2)/(4*D)-0.001;
				#endif 

				/* To satisfy stability conditions, at given dx, dt repeat central approximation finite difference diffusion until reach 30 min tick*/
				for (int i = 0; i < 30/deltat; i++) this->diffuseChem(ichem, deltat, D);
			}
		#endif // GPU_DIFFUSE
	#else // PDE_DIFFUSE
		this->NetlogoDiffuse();
	#endif // PDE_DIFFUSE
}

void WHWorld::runCells() {
	int cellsSize = cells.size(); /* This is only an upper bound on cell list size. It is NOT an actual count of cells (some entries are NULL) */
	#pragma omp parallel for
	for (int i = 0; i < cellsSize; i++) {
		Cell* cell = cells.getDataAt(i);
		if (!cell) continue;
		cell->cellFunction();
	}
}

void WHWorld::executeCells() {
	#ifdef PROFILE_CELL_FUNC
		TIME_STAGE(this->runCells(),		"Cell Function: Chondrocytes",	"	");
	#else
    	cout << " execute cells " << endl;
		this->runCells();
	#endif
}

void WHWorld::executeECMs(){
	cerr << " ECM function  " << endl;
	int numPatches = (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny;
	#pragma omp parallel for
	for (int in = 0; in < numPatches; in++) {
		if (worldECM[in].empty[read_t] == false) this->worldECM[in].ECMFunction();
	}
}

void WHWorld::requestECMfragments() {
	if (WHWorld::highTNFdamage == true) { 
		cout << " high TNF damage " << endl;
		WHWorld::highTNFdamage = false;
		for (int in = 0; in < (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny; in++) {
			#ifndef CALIBRATION
				if (this->WHWorldChem.pTNF[in] > 10) { 
			#else
				if (this->WHWorldChem.pTNF[in] > 10) {
			#endif
					cout << " Degrade ECM " << endl;
					this->worldECM[in].fragmentNCollagen();
					this->worldECM[in].fragmentNAggrecan();
				}
		}
	}
}

/* Each patch diffuses 50% of its chemical equally to its 8 neighboring patches. (Each neighbor receives 1/8 of 50% of the patch's original amount of chemical neighboring patch. */
void WHWorld::NetlogoDiffuse() {
	cerr << " NetLogoDiffuse " << endl;
	for(int ichem = TNF; ichem <= IL1beta; ichem++) {
		for (int iz = 0; iz < nz; iz++) {
			for (int iy = 0; iy < ny; iy++) {
				for (int ix = 0; ix < nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;  // Patch row major index
					if (this->chemAllocation[ichem][in] <= 0) {
						this->chemAllocation[ichem][in] = 0;
						continue;  // No chemical to diffuse
					}

					float value = 0.5*this->chemAllocation[ichem][in]/8;
					for (int dx = -1; dx < 2; dx++) {
						for (int dy = -1; dy < 2; dy++) {
							int dz = 0;
							int neighborindex = (ix + dx) + (iy + dy)*nx + (iz + dz)*ny*nx;
							if (ix + dx < 0 || ix + dx >= nx || iy + dy < 0 || iy + dy >= ny) continue;
							if (dx == 0 && dy == 0 ) continue;
							if (ix%(nx - 1) == 0 && dx == 1) continue;  // TODO(Kim): INSERT REF? (What is this?)
							if (ix%nx == 0 && dx == -1) continue;  // TODO(Kim): INSERT REF? (What is this?)
							this->chemAllocation[ichem + 8][neighborindex] += value;
							this->chemAllocation[ichem + 8][in] -= value;
						}
					}
				}
			}
		}
	}
	for (int xi = 0; xi < nx; xi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int zi = 0; zi < nz; zi++) {
				int in = xi + yi*nx + zi*nx*ny;

				// Update patch chemical concentration:
				this->WHWorldChem.pTNF[in] = this->WHWorldChem.dTNF[in] + this->WHWorldChem.pTNF[in];
				this->WHWorldChem.pTGF[in] = this->WHWorldChem.dTGF[in] + this->WHWorldChem.pTGF[in];
				this->WHWorldChem.pIL1beta[in] = this->WHWorldChem.dIL1beta[in] + this->WHWorldChem.pTGF[in];
				this->WHWorldChem.dTNF[in] = 0;
				this->WHWorldChem.dTGF[in] = 0;
				this->WHWorldChem.dIL1beta[in] = 0;
			}
		}
	}
}

#ifdef GPU_DIFFUSE
	void WHWorld::diffuseChemGPU(){
		cerr << "Diffuse Chem (GPU)" << endl;
		// Loop over all types of chemical and perform convolution-based diffusion on GPU
		const int num_basechem_types = this->baselineChem.size();
		for (int ic = 0; ic < num_basechem_types; ic++) {
			//cerr << "	Diffusing type " << ic << " of " << num_basechem_types - 1 << endl;
			if (!computeChemDiffusionGPU(
				this->h_diffusion_results[ic],   // output: t<chem>
				this->chemAllocation[ic],		 // input:  p<chem>
				this->h_dKernel_spectrum[ic],
				*(this->chem_cctx),
				WHWorld::clock)) {} // TODO: Error handling
		}
	}
#endif	// GPU_DIFFUSE

// Discretization of PDE Diffusion Equation using central difference approximation
// Note: diffuseChem() will very likely get replaced by a new correct version, thus this doesn't need comments just yet
void WHWorld::diffuseChem(int ichem, float dt, float coeff){
	// Calculate change in concentration over dt at each patch
	float* tempPtr = new float[nx*ny*nz]; 
	#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			if (tid == 0) *ntp = omp_get_num_threads();
			start_times[tid] = omp_get_wtime();
			#pragma omp for nowait                                                                                                                                                                                                                                                                                                                                                                  
	#else
		#ifdef V_a
			#pragma omp parallel for schedule(dynamic)
		#else
			#pragma omp parallel
				{
			#pragma omp for
		#endif	// V_a or V_b
	#endif

    // Calculate central difference approximation at each patch along each direction
    for(int yi = 1; yi < ny - 1; yi++){
        for(int xi = 1; xi < nx - 1; xi++){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
			#ifdef MODEL_3D
				for(int zi = 1; zi < nz - 1; zi++){
			#else
                int zi = 0;
			#endif
                    int index = xi+ yi*nx + zi*nx*ny;
                    REAL d;
                    int XPlusOne = (xi+1) + yi*nx + zi*nx*ny;
                    int XMinusOne = (xi-1) + yi*nx + zi*nx*ny;
                    int YPlusOne = xi + (yi+1)*nx + zi*nx*ny;
                    int YMinusOne = xi + (yi-1)*nx + zi*nx*ny;
					#ifdef MODEL_3D
						int ZPlusOne = xi + yi*nx + (zi+1)*nx*ny;
						int ZMinusOne = xi + yi*nx + (zi-1)*nx*ny;
					#endif

                    // Central Difference Approximation along x and y directions:                          
                    float d2phi_dx2 = (this->chemAllocation[ichem][XPlusOne] - 2.*this->chemAllocation[ichem][index] + this->chemAllocation[ichem][XMinusOne])/(dx*dx);
                    float d2phi_dy2 = (this->chemAllocation[ichem][YPlusOne] - 2.*this->chemAllocation[ichem][index] + this->chemAllocation[ichem][YMinusOne])/(dy*dy);
                    // 2D Central Differnce Approximation of Diffusion:
                    tempPtr[index] = dt*coeff*(d2phi_dx2 + d2phi_dy2);
					#ifdef MODEL_3D
						// Central Difference Approximation along z direction                          
						float d2phi_dz2 = (this->chemAllocation[ichem][ZPlusOne] - 2.*this->chemAllocation[ichem][index] + this->chemAllocation[ichem][ZMinusOne])/(dz*dz);
						// 3D Central Differnce Approximation of Diffusion 
						tempPtr[index] = dt*coeff*(d2phi_dx2 + d2phi_dy2 + d2phi_dz2);
				}
					#endif
        }
    }

	#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
		end_times[tid] = omp_get_wtime();
		elapsed1[tid] = end_times[tid] - start_times[tid];
		#pragma omp barrier
	#endif

	#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
		start_times[tid] = omp_get_wtime();
		#pragma omp for nowait
	#else
		#ifdef V_a
			#pragma omp parallel for schedule(dynamic)
		#else
			#pragma omp for
		#endif	// V_a or V_b
	#endif
	
	// Update concentration from central difference scheme
	for(int yi=0; yi<ny; yi++){
		#pragma omp simd
		for(int xi=0; xi<nx; xi++){
			#ifdef MODEL_3D
				for(int zi=0; zi<nz; zi++){
			#else
                int zi = 0;
			#endif
            
			int index = xi+ yi*nx + zi*nx*ny;
            if (yi == 0 || yi == ny-1 || xi == 0 || xi == nx-1  || zi == 0 || zi == nz-1) { // constant padding boundary condition
                int countCaAlg = WHWorld::initialCaAlg;
                this->chemAllocation[ichem][index] = this->baselineChem[ichem]/countCaAlg;
            } else {
                this->chemAllocation[ichem][index] += tempPtr[index];
            }

			#ifdef MODEL_3D
				}
			#endif
		}
	}
        }
		#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
			end_times[tid] = omp_get_wtime();
			elapsed2[tid] = end_times[tid] - start_times[tid];
			#pragma omp barrier
			}

			cout << "Chemical " << ichem << ":" << endl;
			for(int t = 0; t < num_threads; t++) cout << "	thread " << t << " took: [" << elapsed1[t] << "]	[" << elapsed2[t] << "]" << endl;	
		#endif
	delete[] tempPtr;
}

void WHWorld::updateTotalChem(){
	this->WHWorldChem.totalTNF = 0;
	this->WHWorldChem.totalTGF = 0;
	this->WHWorldChem.totalIL1beta = 0;

	float sumTNF = 0, sumTGF = 0, sumIL1 = 0;
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				int in = xi + yi*nx + zi*nx*ny;
				this->WHWorldChem.pTNF[in] = this->WHWorldChem.dTNF[in] + this->WHWorldChem.pTNF[in];
				this->WHWorldChem.pTGF[in] = this->WHWorldChem.dTGF[in] + this->WHWorldChem.pTGF[in];
				this->WHWorldChem.pIL1beta[in] = this->WHWorldChem.dIL1beta[in] + this->WHWorldChem.pIL1beta[in];

				this->WHWorldChem.dTNF[in] = 0;
				this->WHWorldChem.dTGF[in] = 0;
				this->WHWorldChem.dIL1beta[in] = 0;

				// Update gradient
				float patchIL1beta = this->WHWorldChem.pIL1beta[in];
				float patchTNF = this->WHWorldChem.pTNF[in];
				float patchTGF = this->WHWorldChem.pTGF[in];
				this->WHWorldChem.pcellgrad[in] = patchTGF;

				sumTNF += this->WHWorldChem.pTNF[in];
				sumTGF += this->WHWorldChem.pTGF[in];
				sumIL1 += this->WHWorldChem.pIL1beta[in];
			}
		}
	}
    this->WHWorldChem.totalTNF += sumTNF;
	this->WHWorldChem.totalTGF += sumTGF;
	this->WHWorldChem.totalIL1beta += sumIL1;
}

// Always called in non-GPU_DIFFUSE. Called only at beginning otherwise
void WHWorld::updateChemCPU() {
	int totaldam = countPatchType(damage);
	this->WHWorldChem.totalTNF = 0;
	this->WHWorldChem.totalTGF = 0;
	this->WHWorldChem.totalIL1beta = 0;

    float sumTNF = 0, sumTGF = 0, sumIL1 = 0;
	for (int zi = 0; zi < nz; zi++) {
		#pragma omp parallel for reduction(+:sumTNF, sumTGF, sumIL1) 
			for (int yi = 0; yi < ny; yi++) {
		#pragma omp simd
				for (int xi = 0; xi < nx; xi++) {
					int in = xi + yi*nx + zi*nx*ny;
					
					// Update patch chemical concentration
					#ifndef CALIBRATION
						this->WHWorldChem.pTNF[in] = this->WHWorldChem.dTNF[in] + this->WHWorldChem.pTNF[in];//*0.02; //					//this->WHWorldChem.pTNF[in] = this->WHWorldChem.dTNF[in] + (this->WHWorldChem.pTNF[in])*(WHWorld::cytokineDecay[0]);
						this->WHWorldChem.pTGF[in] = this->WHWorldChem.dTGF[in] + this->WHWorldChem.pTGF[in];//*0.02; //					//this->WHWorldChem.pTGF[in] = this->WHWorldChem.dTGF[in] + (this->WHWorldChem.pTGF[in])*(WHWorld::cytokineDecay[1]);
						this->WHWorldChem.pIL1beta[in] = this->WHWorldChem.dIL1beta[in] + this->WHWorldChem.pIL1beta[in];//*0.02; //					//this->WHWorldChem.pIL1beta[in] = this->WHWorldChem.dIL1beta[in] + (this->WHWorldChem.pIL1beta[in])*(WHWorld::cytokineDecay[4]);
					#else
						this->WHWorldChem.pTNF[in] = this->WHWorldChem.dTNF[in] + this->WHWorldChem.pTNF[in]*0.02;
						this->WHWorldChem.pTGF[in] = this->WHWorldChem.dTGF[in] + this->WHWorldChem.pTGF[in]*0.02;
						this->WHWorldChem.pIL1beta[in] = this->WHWorldChem.dIL1beta[in] + this->WHWorldChem.pIL1beta[in]*0.02;
					#endif
					
					this->WHWorldChem.dTNF[in] = 0;
					this->WHWorldChem.dTGF[in] = 0;
					this->WHWorldChem.dIL1beta[in] = 0;

					// Update gradient
					float patchIL1beta = this->WHWorldChem.pIL1beta[in];
					float patchTNF = this->WHWorldChem.pTNF[in];
					float patchTGF = this->WHWorldChem.pTGF[in];
					this->WHWorldChem.pcellgrad[in] = patchTGF;

					// Update total chemical values
					sumTNF += this->WHWorldChem.pTNF[in];
					sumTGF += this->WHWorldChem.pTGF[in];
					sumIL1 += this->WHWorldChem.pIL1beta[in];
				}
			}
	}
    this->WHWorldChem.totalTNF += sumTNF;
	this->WHWorldChem.totalTGF += sumTGF;
	this->WHWorldChem.totalIL1beta += sumIL1;
}

void WHWorld::updateChem() {
	#ifdef GPU_DIFFUSE
		int totaldam = countPatchType(damage);
		this->WHWorldChem.totalTNF = 0;
		this->WHWorldChem.totalTGF = 0;
		this->WHWorldChem.totalIL1beta = 0;

		float sumTNF = 0, sumTGF = 0, sumIL1 = 0;
		for (int zi = 0; zi < nz; zi++) {
			#pragma omp parallel for reduction(+:sumTNF, sumTGF, sumIL1)
				for (int yi = 0; yi < ny; yi++) {
			#pragma omp simd
					for (int xi = 0; xi < nx; xi++) {
						int in = xi + yi*nx + zi*nx*ny;
						
						// Update patch chemical concentration
						this->WHWorldChem.pTNF[in] = this->WHWorldChem.dTNF[in] + this->WHWorldChem.tTNF[in]*0.02;
						this->WHWorldChem.pTGF[in] = this->WHWorldChem.dTGF[in] + this->WHWorldChem.tTGF[in]*0.02;
						this->WHWorldChem.pIL1beta[in] = this->WHWorldChem.dIL1beta[in] + this->WHWorldChem.tIL1beta[in]*0.02;

						this->WHWorldChem.dTNF[in] = 0;
						this->WHWorldChem.dTGF[in] = 0;
						this->WHWorldChem.dIL1beta[in] = 0;

						// Update gradient
						float patchIL1beta = this->WHWorldChem.pIL1beta[in];
						float patchTNF = this->WHWorldChem.pTNF[in];
						float patchTGF = this->WHWorldChem.pTGF[in];
						this->WHWorldChem.pcellgrad[in] = patchTGF;

						// Update total chemical values
						sumTNF += this->WHWorldChem.pTNF[in];
						sumTGF += this->WHWorldChem.pTGF[in];
						sumIL1 += this->WHWorldChem.pIL1beta[in];
					}
				}
		}
			this->WHWorldChem.totalTNF += sumTNF;
			this->WHWorldChem.totalTGF += sumTGF;
			this->WHWorldChem.totalIL1beta += sumIL1;
			cout << "		Total TNF: " << sumTNF << ", Total TGF: " << sumTGF << ", Total IL1beta: " << sumIL1 << endl;
	#else
		// Always calling updateChemCPU() for non-GPU_DIFFUSE chem updates
		updateChemCPU();
	#endif
}

void WHWorld::executeAllECMUpdates() {
	for (int iz = 0; iz < nz; iz++) {
		#pragma omp parallel for
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
				this->worldECM[in].updateECM();
			}
		}
	}
}

void WHWorld::executeAllECMResetRequests() {
	for (int iz = 0; iz < nz; iz++) {
		#pragma omp parallel for
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
				this->worldECM[in].resetrequests();
			}
		}
	}
}

void WHWorld::updateECMManagers() {
	#ifdef PROFILE_ECM_UPDATE
		TIME_STAGE(this->executeAllECMUpdates(),		"	updateECM()",		"	");
		TIME_STAGE(this->executeAllECMResetRequests(),	"	resetrequests()",		"	");
	#else
		this->executeAllECMUpdates();
		this->executeAllECMResetRequests();
	#endif
}

void WHWorld::updatePatches() {
	for (int iz = 0; iz < nz; iz++) {
		#pragma omp parallel for
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
				this->worldPatch[in].updatePatch();
			}
		}
	}
}

/*
 * Steps:
 * 1. Perform updates
 * 2. Remove all dead cells
 * 3. If OMP, add cells from thread-local lists to corresponding global lists
 */
void WHWorld::updateCells() {
	cerr << "	removing dead cells" << endl;
	int cellsSize = cells.size();
	#pragma omp parallel for
		for (int i = 0; i < cellsSize; i++) {
	#ifdef _OMP
		int tid = omp_get_thread_num();
	#else
		int tid = DEFAULT_TID;
	#endif
		// Get pointer of cell i from the array chain
		Cell* cell = cells.getDataAt(i);
		if (!cell) continue;  // cell was deleted
		cell->updateAgent();

		int DeletedCell = 0;
		// Remove dead cells
		if (cell->isAlive() == false) {
			// Get residing patch index and update its occupancy
			int in = cell->getIndex();
			this->worldPatch[in].clearOccupied();
			this->worldPatch[in].occupiedby[write_t] = nothing;
			this->worldPatch[in].dirty = true;
			cells.deleteData(i, tid);
			delete cell;
			DeletedCell++;
			/* Added by MM to check types of cell stages and subtract from respective counters: */
			if (typeid(*this) == typeid(Stem)) {
				Stem::numOfStem--;
			}
			else if (typeid(*this) == typeid(Progen)) {
				Progen::numOfProgen--;
			}
			else if (typeid(*this) == typeid(NP)) {
				NP::numOfNP--;
			}
		}
	}
	Cell::numOfCells = cells.actualSize();

	// Add new cells
	#ifdef _OMP
	/* In OMP, cells were only added to each thread's local list when sproutAgent() was called. 
	* Thus, this step is needed to add those cells onto the global lists. */
		cerr << "	updateCell() _OMP" << endl;
		// TODO: parallelize
		//int numThreads = omp_get_num_threads();
		int numThreads = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
		//cout << "		numThreads = " << numThreads << endl;
		for (int tid = 0; tid < numThreads; tid++) {
			// Cells
			vector<Cell*>* fvec_ptr = localNewCells[tid];
			for (vector<Cell*>::iterator cell_it = fvec_ptr->begin(); cell_it != fvec_ptr->end(); cell_it++) {
				Cell* newCell = *cell_it;
				if (!cells.addData(newCell, tid)) {
					cerr << "Error: Could not add cell" << endl;
					exit(-1);
				}
			}
			fvec_ptr->clear();
		}
	#endif
}

/****************************************************************
 * MAJOR SECTION SUBROUTINES - end                              *
 ****************************************************************/
//NOTE: only use this function to sprout new_coll, new_agg in initialization.

void WHWorld::sproutAgent(int num, int patchType, int agentType,
	int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) {
	#ifdef OPT_CELL_SEEDING
        if (xmin != 0 || xmax != nx || ymin != 0 || ymax != ny || zmin != 0 || zmax != nz) sproutAgentInArea (num, patchType, agentType, xmin, xmax, ymin, ymax, zmin, zmax);
        else if (patchType == CaAlg) sproutAgentInArea (num, patchType, agentType, xmin, xmax, ymin, ymax, zmin, zmax); 
        else sproutAgentInWorld (num, patchType, agentType); 
	#else
  		// Target a specific area of the world
		sproutAgentInArea (num, patchType, agentType, xmin, xmax, ymin, ymax, zmin, zmax);
	#endif
}

void WHWorld::sproutAgentInArea(int num, int patchType, int agentType, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) {
	int count = 0;
	vector <int> patchlist;
	int* reservoir = new int [num];
	for (int i = 0; i < num; i++) reservoir[i] = -1;
	Patch* tempPatchPtr;
	Agent* tempAgentPtr;
	int in, agentIndex, max;
	for (int izz = zmin; izz < zmax + 1; izz++) {
		for (int iyy = ymin; iyy < ymax + 1; iyy++) {
			for (int ixx = xmin; ixx < xmax + 1; ixx++) {
				in = ixx + iyy*nx + izz*nx*ny;
        		
				// Try another patch if this one is out of bounds or the wrong type or occupied
				if (ixx < 0 || ixx >= nx || iyy < 0 || iyy >= ny || izz < 0 || izz >= nz) continue;
				if (WHWorld::worldPatch[in].type[read_t] != patchType) continue;
				if (this->worldPatch[in].isOccupied() == false) patchlist.push_back(in);
			}
		}
	}
	for (int i = 0; i < num; i++) {
		if (patchlist.size() == 0) {  // No available patches
			cout << " sprout agent error, no available patch within bounds! " << endl;
			delete[] reservoir;
			return;
		}
		int randnumber = rand() % patchlist.size();
		reservoir[i] = patchlist[randnumber];  // Prepare 'num' random patches
	}

	// Sprout agent on each patch in reservoir
	for (int i = 0; i < num; i++) {
		int in = reservoir[i];
		if (in < 0 || in > (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny) continue;
			switch (agentType) {
				case stem: {
					cout << "new stem added" << endl; //added for debugging
					tempPatchPtr = &(this->worldPatch[in]);
					Stem* newStem = new Stem(tempPatchPtr);
					#ifdef _OMP
						int tid = omp_get_thread_num();
						this->localNewCells[tid]->push_back(newStem);
					#else
						if (!this->cells.addData(newStem, DEFAULT_TID)) {
							cerr << "Error: Could not add stem cell in sproutAgent()" << endl;
							exit(-1);
						}
					#endif
					/* Added by MM to check types of cell stages and add to respective counters: */
					if (typeid(*this) == typeid(Stem)) {
						Stem::numOfStem++;
					}
					Cell::numOfCells++;
        			
					this->worldPatch[in].setOccupied();
        			this->worldPatch[in].occupiedby[write_t] = stem;
					this->worldPatch[in].dirty = true;
					cout << "patch index " << in << endl; //added for debugging
        			
					break;
      			}
				case progen: {
					cout << "new progen added" << endl; //added for debugging
					tempPatchPtr = &(this->worldPatch[in]);
					Progen* newProgen = new Progen(tempPatchPtr);
					#ifdef _OMP
						int tid = omp_get_thread_num();
						this->localNewCells[tid]->push_back(newProgen);
					#else
						if (!this->cells.addData(newProgen, DEFAULT_TID)) {
							cerr << "Error: Could not add pre-np cell in sproutAgent()" << endl;
							exit(-1);
						}
					#endif
					/* Added by MM to check types of cell stages and add to respective counters: */
					if (typeid(*this) == typeid(Progen)) {
						Progen::numOfProgen++;
					}
					Cell::numOfCells++;

					this->worldPatch[in].setOccupied();
					this->worldPatch[in].occupiedby[write_t] = progen;
					this->worldPatch[in].dirty = true;
					
					break;
				}
				case np: {
					cout << "new np added" << endl; //added for debugging
					tempPatchPtr = &(this->worldPatch[in]);
					NP* newNP = new NP(tempPatchPtr);
					#ifdef _OMP
						int tid = omp_get_thread_num();
						this->localNewCells[tid]->push_back(newNP);
					#else
						if (!this->cells.addData(newNP, DEFAULT_TID)) {
							cerr << "Error: Could not add np cell in sproutAgent()" << endl;
							exit(-1);
						}
					#endif
					/* Added by MM to check types of cell stages and add to respective counters: */
					if (typeid(*this) == typeid(NP)) {
						NP::numOfNP++;
					}
					Cell::numOfCells++;

					this->worldPatch[in].setOccupied();
					this->worldPatch[in].occupiedby[write_t] = np;
					this->worldPatch[in].dirty = true;

					break;
				}
      			case orig_coll: {
					this->worldECM[in].ocollagen[write_t] = this->worldECM[in].ocollagen[read_t] + 1;
					this->worldECM[in].dirty = true;
					this->worldECM[in].isEmpty();
					break;
				}
				case orig_agg: {
					this->worldECM[in].oaggrecan[write_t] = this->worldECM[in].oaggrecan[read_t] + 1;
					this->worldECM[in].dirty = true;
					this->worldECM[in].isEmpty();
					break;
				}
				case new_coll: {
					this->worldPatch[in].initcollagen = true;
					this->worldECM[in].ncollagen[write_t] = this->worldECM[in].ncollagen[read_t] + 1;
					this->worldECM[in].dirty = true;
					this->worldECM[in].isEmpty();
					break;
				}
				case new_agg: {
					this->worldPatch[in].initaggrecan = true;
					this->worldECM[in].naggrecan[write_t] = this->worldECM[in].naggrecan[read_t] + 1;
					this->worldECM[in].dirty = true;
					this->worldECM[in].isEmpty();
					break;
				}
			}
		}
	delete[] reservoir;
	return;
}

#ifdef OPT_CELL_SEEDING
	/*
	* Optimized by:
	*  - If sprout in tissue:
	*     (*) Randomly choosing a target patch:
	*          - If is tissue and occupied, sprout
	* 			    - Else repeat (*)
	*  - Else (sprout in blood):
	*     (**)	Look at the list of capillary patches initialized in the setup stage
	*          - Create a list of unoccupied capillary patches
	*          - Pick randomly and sprout
	*          - Repeat until 'num' cells are sprouted
	*/
	void WHWorld::sproutAgentInWorld(int num, int patchType, int agentType) {// bool bloodORtiss
		int count = 0;
		vector <int> patchlist;
		int* reservoir = new int [num];
		for (int i = 0; i < num; i++) reservoir[i] = -1;
		Patch* tempPatchPtr;
		Agent* tempAgentPtr;
		int in, agentIndex, max;
		int totalNumPatches = this->(n-1)x*this->(ny-1)*this->(nz-1);	 // Is this accurate?
		int numfound;
		int counter;
		int threshold;
		vector<int> unoccupiedCaps;
		switch (patchType) {
			for (int i = 0; i < num; i++) {
				if (reservoir[i] < 0 || reservoir[i] > (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny) continue;
				switch (agentType) {
					case stem: {
						tempPatchPtr = &(this->worldPatch[reservoir[i]]);
						Stem* newStem = new Stem(tempPatchPtr);
						#ifdef _OMP
							int tid = omp_get_thread_num();
							this->localNewCells[tid]->push_back(newStem);
						#else
							if (!this->cells.addData(newStem, DEFAULT_TID)) {
								cerr << "Error: Could not add stem cell in sproutAgentInWorld()" << endl;
								exit(-1);
							}
						#endif
						//this->worldPatch[reservoir[i]].occupied = true;
						this->worldPatch[reservoir[i]].setOccupied();
						this->worldPatch[reservoir[i]].occupiedby = stem;
						break;
					}
					case progen: {
						tempPatchPtr = &(this->worldPatch[reservoir[i]]);
						Progen* newProgen = new Progen(tempPatchPtr);
						#ifdef _OMP
							int tid = omp_get_thread_num();
							this->localNewCells[tid]->push_back(newProgen);
						#else
							if (!this->cells.addData(newProgen, DEFAULT_TID)) {
								cerr << "Error: Could not add pre-np cell in sproutAgentInWorld()" << endl;
								exit(-1);
							}
						#endif
						//this->worldPatch[reservoir[i]].occupied = true;
						this->worldPatch[reservoir[i]].setOccupied();
						this->worldPatch[reservoir[i]].occupiedby = progen;
						break;
					}
					case np: {
						tempPatchPtr = &(this->worldPatch[reservoir[i]]);
						NP* newNP = new NP(tempPatchPtr);
						#ifdef _OMP
							int tid = omp_get_thread_num();
							this->localNewCells[tid]->push_back(newNP);
						#else
							if (!this->cells.addData(newNP, DEFAULT_TID)) {
								cerr << "Error: Could not add np cell in sproutAgentInWorld()" << endl;
								exit(-1);
							}
						#endif
						//this->worldPatch[reservoir[i]].occupied = true;
						this->worldPatch[reservoir[i]].setOccupied();
						this->worldPatch[reservoir[i]].occupiedby = np;
						break;
					}
				}
			}
			delete[] reservoir;
			return;
		}
#endif  //OPT_CELL_SEEDING

int WHWorld::countPatchType(int whichType) {
	if (whichType == CaAlg) {
		Patch::numOfEachTypes[whichType] = 0;
		for (int iz = 0; iz < this->nz; iz++) {
			int currCount = 0;
			#pragma omp parallel for reduction(+:currCount)

			for (int iy = 0; iy < this->ny; iy++) {
				for (int ix = 0; ix < this->nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;
					if (this->worldPatch[in].type[read_t] == whichType) {
						currCount++;
						//Patch::numOfEachTypes [whichType]++;
					}
				}
			}
			Patch::numOfEachTypes [whichType] += currCount;
		}
	} else if (whichType == damage) {
		Patch::numOfEachTypes[whichType] = 0;
		for (int iz = 0; iz < this->nz; iz++) {
			int currCount = 0;
			#pragma omp parallel for reduction(+:currCount)
			for (int iy = 0; iy < this->ny; iy++) {
				for (int ix = 0; ix < this->nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;
					currCount += this->worldPatch[in].damage[read_t];
				}
			}
			Patch::numOfEachTypes [whichType] += currCount;
		}
	} else
		cout << "type must be 0, 1, 2 , 3 or 4!" << endl;	//cout << "type must be 0, 1, 2 , 3, 4 or 5!" << endl;
	return Patch::numOfEachTypes[whichType];
}

int WHWorld::mmToPatch(double mm) {
	return mm*(this->patchpermm);
}

int WHWorld::reportTick(int hour, int day) {
	return (hour*2 + day*48);
}

double WHWorld::reportMinute() {
	return (WHWorld::clock)*30;
}

double WHWorld::reportHour() {
	return (WHWorld::clock)/2;
}

double WHWorld::reportDay() {
	return (WHWorld::clock)/48;
}

int WHWorld::countNeighborPatchType(int ix, int iy, int iz,  int patchType) {
	int neighborcount = 0;
	for (int dx = -1; dx < 2; dx++) {
		for (int dy = -1; dy < 2; dy++) {
			for (int dz = -1; dz < 2; dz++) {
				int neighborindex = (ix + dx) + (iy + dy)*nx + (iz + dz)*ny*nx;
				if (ix + dx < 0 || ix + dx >= nx || iy + dy < 0 || iy + dy >= ny || iz + dz < 0 || iz + dz >= nz) continue;
				if (dx == 0 && dy ==0 && dz ==0 ) continue;
				if (Agent::agentPatchPtr[neighborindex].type[read_t] == patchType) neighborcount++;
			}
		}
	}
	return neighborcount;
}

#ifdef MODEL_SCAFFOLD
	void WHWorld::updateSwellingRatio(){
		float Alg_ww = this->Alg_wv/(this->Alg_wv);
		double tmin = reportMinute(); 

		/* Calculate Swelling Ratio (Q) given Alg concentration (% w/w) at a given time (in minutes)
		* 		Q = (a * Alg_ww +  b) * t_min + (c * Alg_ww + d)
		*
		*       Swelling ratio favorable for cell adhesion, growth, diffusion of nutrients
		*       Depends on Alg content (% w/w) 
		*       Rapid swelling in initial 10 min, with slight increase until ~24h
		*/
		#ifdef CALIBRATION
			this->Q = (this->SwellRatio[0]*Alg_ww + this->SwellRatio[1])*log(tmin) + (this->SwellRatio[2]*Alg_ww + this->SwellRatio[3]); 
		#else
			this->Q = (0.4*Alg_ww + 0.4)*log(tmin) + (3*Alg_ww + 7.9); 
		#endif
		this->Q=0;
		//cout << " this->Q =" << (this->SwellRatio[0]<<"*"<<Alg_ww<< " + "<< this->SwellRatio[1])<<"*log("<<tmin<<") + ("<<this->SwellRatio[2]<<"*"<<Alg_ww<< " + "<<this->SwellRatio[3]<<")" << endl; 
		//cout << " Swelling Ratio: " << this->Q << endl; 
	}

	void WHWorld::updateMassLoss(){
		float Alg_ww = this->Alg_wv/(this->Alg_wv);
		float tweek = reportDay()/7;
		float w_t;

		/* Calculate Mass Loss (w_t) of gel with given Alg concentration (% w/w) at current time (in weeks)
		* 		w_t = (a * Alg_ww +  b) * t_weeks + (c * Alg_ww + d)
		*
		*       Degree of dregradation depends greatly on Alg (% w/w) content 
		*       Highest weight loss percentage occurs during first week in vitro, with little weight loss over next 3 weeks 
		*/
		#ifdef CALIBRATION
			w_t = this->MassLoss[0] + this->MassLoss[1]*(pXL) + this->MassLoss[2]*(reportDay()) - this->MassLoss[3]*(pXL)*(reportDay());
		#else
			w_t = 0.234 + 7.785*(pXL) + 0.15*(reportDay()) - 1.36*(pXL)*(reportDay()); //(17.6*Alg_ww - 0.9)*log(tweek) + (60*Alg_ww + 5.3);
		#endif
		
		if (w_t < 0 ) w_t = 0;  // no negative mass loss

		// If there is % mass loss since last call, "degrade" % CaAlg patches and replace with tissue
		if (w_t > this->w){
			float changeInPatches = (0.01)*(w_t - this->w)*WHWorld::initialCaAlg; 
			this->degradeCaAlg(changeInPatches); 
		}
		
		this->w = w_t; 
		//cout << " w_t = "<<this->MassLoss[0]<<" + "<< this->MassLoss[1]<<"*"<<(pXL)<<" + "<<this->MassLoss[2]<<"*"<<(reportDay()) <<" - " <<this->MassLoss[3]<<"*"<<(pXL)<<"*"<<(reportDay()) << endl; 
		//cout << " Mass Loss (%): " << this->w << endl; 
		//cout << " Number of Ca-Alg patches: " << this->countPatchType(CaAlg) << endl;
	}
#endif //MODEL_SCAFFOLD

#ifdef MODEL_SCAFFOLD
void WHWorld::degradeCaAlg(int numOfPatches){
	int xmin = 0; 
	int xmax = nx;
	int ymin = 0; 
	int ymax = ny;
	int zmin = 0; 
	int zmax = nz; 
	int patchType = CaAlg; 
	vector <int> patchlist;
	int* reservoir = new int [numOfPatches];
	for (int i = 0; i < numOfPatches; i++) reservoir[i] = -1;
	int in, agentIndex, max;
    int count = 0; 

	// Make list of possible CaAlg Patches to degrade
	for (int iz = zmin; iz < zmax ; iz++) {
		for (int iy = ymin; iy < ymax ; iy++) {
			for (int ix = xmin; ix < xmax ; ix++) {
				in = ix + iy*nx + iz*nx*ny;
        		
				// Try another patch if this one is out of bounds or the wrong type or occupied
				if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) continue;
				if (WHWorld::worldPatch[in].type[read_t] != patchType) continue;
				patchlist.push_back(in);
			}
		}
	}

    // Choose random patches from patch list 
	for (int i = 0; i < numOfPatches; i++) {
		if (patchlist.size() == 0) {  // No available patches
			//cout << " CaAlg degrade error, no available patch within bounds! " << endl;
			delete[] reservoir;
			return;
		}
		int randnumber = rand() % patchlist.size();
		reservoir[i] = patchlist[randnumber];  // Prepare 'num' random patches
	}        

	// Degrade 'numOfPatches" number of CaAlg patches in reservoir list
	for (int i = 0; i < numOfPatches; i++) {
		int in = reservoir[i];
		if (in < 0 || in > (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny) continue;
		WHWorld::worldPatch[in].type[write_t] = nothing; 
		WHWorld::worldPatch[in].color[write_t] = cnothing;
        WHWorld::worldPatch[in].dirty = true;
        count ++; 
	}
	delete[] reservoir;
}
#endif //MODEL_SCAFFOLD

void WHWorld::debugInfo() {
	int stemSize = 0; int progenSize = 0; int npSize = 0;
	int cellsSize = cells.size();
	for (int i = 0; i < cellsSize; i++) {
		Cell* cell = cells.getDataAt(i);
		if (!cell) continue;
		if (cell->isAlive() == false) continue;
		//if (cell->activate[read_t] == false) f++;
		if (typeid(cell) == typeid(Stem)) {
			stemSize++;
		}
		else if (typeid(cell) == typeid(Progen)) {
			progenSize++;
		}
		else if (typeid(cell) == typeid(NP)) {
			npSize++;
		}
		//else af++;
	}

	int numCaAlg = 0;
	numCaAlg = countPatchType(CaAlg);
	cout << " total patches: " << numCaAlg << endl;
	cout << " total cells: " << cells.actualSize() << endl;
	cout << " stem cells: " << stemSize << endl;
	cout << " pre-np cells: " << progenSize << endl;
	cout << " np cells: " << npSize << endl;
}

/*
 * Steps:
 *   1. If OMP, add cells from thread-local lists to corresponding global lists
 *   2. Perform updates
 */
void WHWorld::updateCellsInitial() {
	// Cell lists should be empty
	// Add new cells
	#ifdef _OMP
		cerr << "	updateCell() _OMP" << endl;
		// TODO: parallelize
		
		int numThreads = omp_get_num_threads();
		for (int tid = 0; tid < numThreads; tid++) {
			/* ------------------------------ Cells ------------------------------ */
			vector<Cell*>* fvec_ptr = localNewCells[tid];
			for (vector<Cell*>::iterator cell_it = fvec_ptr->begin(); cell_it != fvec_ptr->end(); cell_it++) {
				Cell* newCell = *cell_it;
				if(!cells.addData(newCell, tid)) {
					cerr << "Error: Could not add cell" << endl;
					exit(-1);
				}
			}
			fvec_ptr->clear();
		}
	#endif

	/* ------------------------------ Cells ------------------------------ */
	// No need for deletion since these are new cells
	int cellsSize = cells.size();
	#pragma omp parallel for
	for (int i = 0 ; i < cellsSize; i++) {
		#ifdef _OMP
			int tid = omp_get_thread_num();
		#else
			int tid = DEFAULT_TID;
		#endif
		Cell* cell = cells.getDataAt(i);
		if (!cell) continue;
		cell->updateAgent();
		/* Added by MM to check types of cell stages and add to respective counters: */
		if (typeid(*this) == typeid(Stem)) {
			Stem::numOfStem++;
		}
		else if (typeid(*this) == typeid(Progen)) {
			Progen::numOfProgen++;
		}
		else if (typeid(*this) == typeid(NP)) {
			NP::numOfNP++;
		}
	}
	Cell::numOfCells = cells.actualSize();
}

int WHWorld::userInput() {
	// Read input parameters from user-specified file
	ifstream infile(util::getInputFileName());

	int numChem = -1;
	// TODO: Make this check for specific tag (field name)

	if (infile.is_open()) {
		char garbage[100];
        /* -------------------------------- CHEMICALS ------------------------------- */
		//cout << "Reading the number of baseline chemicals..." << endl;
		float temp;
		infile >> garbage;
		infile >> temp;
		//cout << garbage;
		//cin >> garbage;

		numChem = temp;
		this->baselineChem.resize(temp);
		cout << "The number of baseline chemicals are: " << this->baselineChem.size() << endl;
		for (int ichem = 0; ichem < baselineChem.size(); ichem++) {
			//cout << "Reading the baselineChemical " << ichem << endl;
			infile >> garbage;
			infile >> this->baselineChem[ichem];
			//cout << garbage;
			//cin >> garbage;
			cout << "baselineChem " << ichem << " is " << this->baselineChem[ichem] << endl;
		}

        /* ---------------------------------- CELLS --------------------------------- */
		//cout << "Reading the initial number of types of cells..." << endl;
		int tempCells;
		infile >> garbage;
		infile >> tempCells;
		//cout << garbage;
		//cin >> garbage;
		this->initialCells.resize(tempCells);
		cout << "The number of types of cells is: " << this->initialCells.size() << endl;
		for (int icell = 0; icell< initialCells.size(); icell++) {
			//cout << "Reading the initial # of cell type " << icell + 1 << endl;
			infile >> garbage;
			infile >> this->initialCells[icell];
			//cout << garbage;
			//cin >> garbage;
			cout<<"# of cell " << icell << " is " << this->initialCells[icell] << endl;
		}

        /* -------------------------- Ca-Alg PROPERTIES -------------------------- */
		//cout << "Reading 1% (w/v) Alg volume (mL)" << endl;
		infile >> garbage;
		infile >> this->Alg_v;
		cout << "Volume of Alg (mL) = " << this->Alg_v << endl;

		//cout << "Reading 1.67% (w/v) Ca volume (mL)" << endl;
		infile >> garbage;
		infile >> this->Ca_v;
		cout << "Volume of Ca (mL) = " << this->Ca_v << endl;

        float totalVolume = this->Alg_v + this->Ca_v; 

        this->Alg_wv = 2;
		//this->Alg_wv = 1.95;

		cout << "Total Volume from file (mL): " << totalVolume << endl;

        /* --------------------------- CYTOKINE PROPERTIES -------------------------- */
		//cout << "Reading Cytokine Properties..." << endl;
		float D;
		int HL_s;
		char tag[200];
		do {
			if (numChem < 1) {
				cerr << "Warning: No chem allocated!!!" << endl;
				break;
			}
			// Allocate memory for diffusion coefficients and half-lifes
			this->D = (float*) malloc(sizeof(float) * numChem);
			this->HalfLifes = (int *) malloc(sizeof(int) * numChem);
			while (infile >> tag) {
				//infile >> tag;
				if (!strcmp(tag, "D:")) {
					infile >> D;
					for (int ic = 0; ic < numChem; ic++) this->D[ic] = D;	
					cout << "	setting all D to: " << D << endl;
				} else if (!strcmp(tag, "HL:")) {
					infile >> HL_s;
					for (int ic = 0; ic < numChem; ic++) this->HalfLifes[ic] = HL_s;	
					cout << "	setting all HLs to: " << HL_s << endl;
				
				} else if (!strcmp(tag, "D_TNF:")) {
					infile >> D;
					this->D[TNF] = D;
					cout << "	D_TNF: " << D << endl;
				} else if (!strcmp(tag, "HL_TNF:")) {
					infile >> HL_s;
					this->HalfLifes[TNF] = HL_s;
					cout << "	HL_TNF: " << HL_s << endl;
				} else if (!strcmp(tag, "D_TGF:")) {
					infile >> D;
					this->D[TGF] = D;
					cout << "	D_TGF: " << D << endl;
				} else if (!strcmp(tag, "HL_TGF:")) {
					infile >> HL_s;
					this->HalfLifes[TGF] = HL_s;
					cout << "	HL_TGF: " << HL_s << endl;
				} else if (!strcmp(tag, "D_IL1beta:")) {
				 	infile >> D;
				 	this->D[IL1beta] = D;
				 	cout << "	D_IL1beta: " << D << endl;
				 } else if (!strcmp(tag, "HL_IL1beta:")) {
				 	infile >> HL_s;
				 	this->HalfLifes[IL1beta] = HL_s;
				 	cout << "	HL_IL1beta: " << HL_s << endl;
				} else  cout << "	invalid tag: " << tag << endl;
			}
		} while (0);

		infile.close();
        cout << "-------------------------------------------" << endl; 
	}  // end of if file opens properly
	else cout << "Cannot open file!" << endl;


	//this->baselineChem[3]; //added manually for testing
	//cout << "The number of baseline chemicals are: " << this->baselineChem.size() << endl;

	//this->initialCells[3]; //added manually for testing
	//cout << "The number of types of cells is: " << this->initialCells.size() << endl;

	//this->Alg_v = 0.0275; //added manually for testing
	//this->Ca_v = 0.005;

	//float totalVolume = this->Alg_v + this->Ca_v;

	//cout << "Total Volume (mL): " << totalVolume << endl;

	return 0;
}

void WHWorld::outputWorld_csv() {
	if (this->clock == 0) {
		//remove( "output/Output_Biomarkers_90Mw_34mM.csv");
		//remove( "output/Output_Biomarkers_90Mw_22mM.csv");
		//remove( "output/Output_Biomarkers_1500Mw_29mM.csv");
        remove( "output/Output_Biomarkers_1500Mw_14mM.csv");
		//remove( "output/Output_Biomarkers_Experiment_200Mw_20mM.csv");
        
		
		//ofstream output_file("output/Output_Biomarkers_Experiment_200Mw_20mM.csv", ios::app);
		ofstream output_file("output/Output_Biomarkers.csv", ios::app);	
		//ofstream output_file("output/Output_Biomarkers_1500Mw_29mM.csv", ios::app);		
		//ofstream output_file("output/Output_Biomarkers_90Mw_22mM.csv", ios::app);		
		//ofstream output_file("output/Output_Biomarkers_90Mw_34mM.csv", ios::app);		

		output_file << "clock (30 min)" << "," << "Day" << "," << "Total TNF (pg)" <<  "," << "Total IL1b (pg)" << "," << "Total TGF (pg)" << "," << "Collagen (ug)" << "," << "Aggrecan (ug)" << "," << "Activated Chondrocytes" << "," << "Total Chondrocytes" << ", Elastic Modulus (kPa) " << ", Swelling Ratio " << ", Mass Loss (%) " << ", Alginate_wv (%)" << ", Alginate_Mw (kDa)" << ", Ca_XL (M)"<< ", Viability Rate (%)" << endl; //output_file << "Tropocollagen" << "," << "Collagen" << "," << "FragentedCollagen" << "," << "Tropoaggrecan" << "," << "Aggrecan" << "," << "FragmentedAggrecan" << "," << "HA" << "," << "FragmentedHA" << "," << "Damage" endl;
        output_file.close();
	}

    //ofstream output_file("output/Output_Biomarkers_Experiment_200Mw_20mM.csv", ios::app);
	ofstream output_file("output/Output_Biomarkers.csv", ios::app);
	//ofstream output_file("output/Output_Biomarkers_1500Mw_29mM.csv", ios::app);
	//ofstream output_file("output/Output_Biomarkers_90Mw_22mM.csv", ios::app);
	//ofstream output_file("output/Output_Biomarkers_90Mw_34mM.csv", ios::app);

	int f = 0; int af = 0;
	int orig_coll = 0; int frag_coll = 0; double new_coll = 0; 
	int orig_agg = 0; int frag_agg = 0; double new_agg = 0; 
	int HA = 0; int fHA = 0;

	int cellsSize = cells.size();
	for (int i = 0; i < cellsSize; i++) {
		Cell* cell = cells.getDataAt(i);
		if (!cell) continue;
		if (cell->isAlive() == false) continue;
		if (cell->activate[read_t] == false) f++;
		if (cell->type[read_t] == stem) {
			Stem::numOfStem++;
		}
		else if (cell->type[read_t] == progen) {
			Progen::numOfProgen;
		}
		else if (cell->type[read_t] == np) {
			NP::numOfNP;
		}
		//if (typeid(cell) == typeid(Stem)) {
		//	stemSize++;
		//}
		//else if (typeid(cell) == typeid(Progen)) {
		//	progenSize++;
		//}
		//else if (typeid(cell) == typeid(NP)) {
		//	npSize++;
		//}
		else af++;
	}

	cout << " total cells: " << cells.actualSize() << endl;
	cout << " stem cells: " << Stem::numOfStem << endl;
	cout << " pre-np cells: " << Progen::numOfProgen << endl;
	cout << " np cells: " << NP::numOfNP << endl;

	for (int in = 0; in < (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny; in++) {
		orig_coll += this->worldECM[in].ocollagen[read_t];
		new_coll += this->worldECM[in].ncollagen[read_t];
		frag_coll += this->worldECM[in].fcollagen[read_t];
		orig_agg += this->worldECM[in].oaggrecan[read_t];
		new_agg += this->worldECM[in].naggrecan[read_t];
		frag_agg += this->worldECM[in].faggrecan[read_t];
	}
	this->countPatchType(damage);
	output_file << this->clock << ",";
	output_file << (this->clock)/48 << ",";
	output_file << this->WHWorldChem.totalTNF << ",";
	output_file << this->WHWorldChem.totalIL1beta << ",";
	output_file << this->WHWorldChem.totalTGF << ",";
	output_file << fixed << std::setprecision(5) << new_coll/1000.0 << "," << new_agg/1000.0 << ","; //ECM 	//output_file << orig_coll << "," << new_coll << "," << frag_coll << "," << orig_agg << "," ; 	//output_file << new_agg << "," << frag_agg << "," << HA << "," << fHA << "," << Patch::numOfEachTypes[4] << "," ;
	//output_file << fixed << std::setprecision(5) << new_coll << "," << new_agg << ","; //ECM 	//output_file << orig_coll << "," << new_coll << "," << frag_coll << "," << orig_agg << "," ; 	//output_file << new_agg << "," << frag_agg << "," << HA << "," << fHA << "," << Patch::numOfEachTypes[4] << "," ;

	//output_file << af << "," << f+af << ","; //cells
	output_file << Stem::numOfStem << "," << Progen::numOfProgen << "," << NP::numOfNP << "," << cells.actualSize() <<","; // cell counts

	#ifdef MODEL_SCAFFOLD
		output_file << this->E << " , " << this->Q << ", " << this->w << "," << this->Alg_wv << ","<< this->Alg_Mn << "," << this->pXL << endl;
	#else
		output_file  << endl;
	#endif
	output_file.close();

	cout << " Collagen: " << new_coll << endl;
	cout << " Aggrecan: " << new_agg << endl;

}

void WHWorld::patchassign_csv() {
	int in = 0;
	int Number = 0;
	for (int iz = 0; iz < nz; iz++) {
		char patchassign[50] = "output/patchassign";
		char cells[50] = "output/cells_read";
		char cells_w[50] = "output/cells_write";
		char initcollagen[50] = "output/initcollagen";
		char initaggrecan[50] = "output/initaggrecan";
		char initHA[50] = "output/initHA";
		char damagezone[50] = "output/damagezone";
		char initialdamage[50] = "output/initialdamage";
		char extension[10] = ".csv";
		char tempNumber[20] = "";
		sprintf (tempNumber, "_t%3.0f_z%d", this->clock, Number);
		strcat(patchassign, tempNumber);
		strcat(patchassign, extension);
		strcat(cells, tempNumber);
		strcat(cells, extension);
		strcat(cells_w, tempNumber);
		strcat(cells_w, extension);
		strcat(initHA, tempNumber);
		strcat(initHA, extension);
		strcat(initcollagen, tempNumber);
		strcat(initcollagen, extension);
		strcat(initaggrecan, tempNumber);
		strcat(initaggrecan, extension);
		strcat(damagezone, tempNumber);
		strcat(damagezone, extension);
		strcat(initialdamage, tempNumber);
		strcat(initialdamage, extension);

		// Patch Assign
		ofstream output_file(patchassign, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldPatch[in].type[read_t] == damage || worldPatch[in].damage[read_t] != 0) {
					output_file << "x";
					continue;
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file << "?";
				}
			}
			output_file << endl;
		}
		output_file.close();

		// initHA
		ofstream output_file1(initHA, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldPatch[in].initHA == true) {
					output_file1 << "u";
					continue;
				}
				if (worldPatch[in].type[read_t] == damage || worldPatch[in].damage[read_t] != 0) {
					output_file1 << "x";
					continue;
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file1 << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file1 << "?";
				}
			}
			output_file1 << endl;
		}
		output_file1.close();

		// Damage Zone
		ofstream output_file2(damagezone, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;

				if (worldPatch[in].inDamzone == true) {
					output_file2 << "z";
					continue;
				}
				if (worldPatch[in].type[read_t] == damage || worldPatch[in].damage[read_t] != 0) {
					output_file2 << "x";
					continue;
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file2 << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file2 << "?";
				}
			}
			output_file2 << endl;
		}
		output_file2.close();

		// Initial Damage
		ofstream output_file3(initialdamage, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldECM[in].oaggrecan[read_t] !=0 && worldPatch[in].damage[read_t] != 0) {
					output_file3 << "g";
					continue;
				}
				if (worldECM[in].oaggrecan[read_t] !=0) {
					output_file3 << "m";
					continue;
				}
				if (worldPatch[in].type[read_t] ==damage || worldPatch[in].damage[read_t] != 0) {
					output_file3 << "x";
					continue;
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file3 << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file3 << "?";
				}
			}
			output_file3 << endl;
		}
		output_file3.close();

		// initCollagen
		ofstream output_file4(initcollagen, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldECM[in].ocollagen[read_t] != 0) {
					output_file4 << "k";
					continue;
				}
				if (worldECM[in].fcollagen[read_t] != 0) {
					output_file4 << "f";
					continue;
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file4 << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file4 << "?";
				}
			}
			output_file4 << endl;
		}
		output_file4.close();

		// initAggrecan
		ofstream output_file5(initaggrecan, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldECM[in].oaggrecan[read_t] !=0) {
					output_file5 << "m";
					continue;
				}
				if (worldECM[in].faggrecan[read_t] !=0) {
					output_file5 << "f";
					continue;
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file5 << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file5 << "?";
				}
			}
			output_file5 << endl;
		}
		output_file5.close();

		// Cells
		ofstream output_file6(cells, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldPatch[in].isOccupied()) {
					if (worldPatch[in].occupiedby[read_t] == stem) {
						output_file6 << "f";
						continue;
					}
					else if (worldPatch[in].occupiedby[read_t] == progen) {
						output_file6 << "g";
						continue;
					}
					else if (worldPatch[in].occupiedby[read_t] == np) {
						output_file6 << "h";
						continue;
					}
				}
				if (worldPatch[in].type[read_t] == nothing) {
					output_file6 << "-";
				}
				if (worldPatch[in].type[read_t] == unidentifiable) {
					output_file6 << "?";
				}
			}
			output_file6 << endl;
		}
		output_file6.close();

		ofstream output_file7(cells_w, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (worldPatch[in].isOccupiedWrite()) {
					if (worldPatch[in].occupiedby[read_t] == stem) {
						output_file6 << "f";
						continue;
					}
					else if (worldPatch[in].occupiedby[read_t] == progen) {
						output_file6 << "g";
						continue;
					}
					else if (worldPatch[in].occupiedby[read_t] == np) {
						output_file6 << "h";
						continue;
					}
				}
				if (worldPatch[in].type[write_t] == nothing) {
					output_file7 << "-";
				}
				if (worldPatch[in].type[write_t] == unidentifiable) {
					output_file7 << "?";
				}
			}
			output_file7 << endl;
		}
		output_file7.close();

		Number++;
	}
}