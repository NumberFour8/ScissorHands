#include "kernel.h"

// Globální proměnné
__device__ __constant__ unsigned int d_invN;
__device__ unsigned int* d_N;
__device__ unsigned int* d_3N;

#ifndef USE_TWISTED
	#include "edwards.h"
#else 
	#include "twisted.h"
#endif

cudaError_t compute(const Aux h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF& coeff,const unsigned int WS)
{
	const int WINDOW_SZ	 = WS;							// Velikost okna
	const int PRECOMP_SZ = (1 << WINDOW_SZ)-1;			// Počet bodů, které je nutné předpočítat
	
	cudaEvent_t start,stop;
	float totalTime = 0;
	void *swQw = NULL,*swPc = NULL,*swAx = NULL;
	gpuErrchk(cudaSetDevice(0));
	
	gpuErrchk(cudaEventCreate(&start));
	gpuErrchk(cudaEventCreate(&stop));

	// Alokace potřebných dat
	cuda_Malloc((void**)&swPc,NUM_CURVES*PRECOMP_SZ*4*MAX_BYTES); // Předpočítané body
	cuda_Malloc((void**)&swQw,NUM_CURVES*4*MAX_BYTES);			  // Pomocný bod
	cuda_Malloc((void**)&swAx,sizeof(Aux));						  // Pomocná struktura
	
	// Pomocná struktura
	cuda_Memcpy(swAx,(void*)&h_input,sizeof(Aux),cudaMemcpyHostToDevice);
	
	// Počáteční body
	VOL digit_t* iter = (digit_t*)swPc;
	for (int i = 0;i < NUM_CURVES;i++){
	   cuda_Memcpy((void*)(iter+0*NB_DIGITS),(void*)initPoints[i].X,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+1*NB_DIGITS),(void*)initPoints[i].Y,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+2*NB_DIGITS),(void*)initPoints[i].Z,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+3*NB_DIGITS),(void*)initPoints[i].T,MAX_BYTES,cudaMemcpyHostToDevice);
	   iter += 4*NB_DIGITS;	
	}

	// Další předpočítané body
	dim3 threadsPerBlock(NB_DIGITS,CURVES_PER_BLOCK);
	
	START_MEASURE(start);
	curvesDbl<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)swQw,(void*)swPc,(void*)swAx);
	for (int i = 1; i < PRECOMP_SZ;++i){ // Tady už je iter nastavené na pozici prvních lichých mocnin
		curvesAdd<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)iter,(void*)swQw,(void*)swPc,(void*)swAx); 
		curvesAdd<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)swQw,(void*)iter,(void*)swPc,(void*)swAx);
		iter += NUM_CURVES*4*NB_DIGITS;
	} 
	STOP_MEASURE("Precomputation phase",start,stop,totalTime);


	// Do swQw nakopírovat neutrální prvek
	iter = (digit_t*)swQw;
	for (int i = 0;i < NUM_CURVES;++i){
		cuda_Memcpy((void*)(iter+0*NB_DIGITS),(void*)neutral->X,MAX_BYTES,cudaMemcpyHostToDevice);
		cuda_Memcpy((void*)(iter+1*NB_DIGITS),(void*)neutral->Y,MAX_BYTES,cudaMemcpyHostToDevice);
		cuda_Memcpy((void*)(iter+2*NB_DIGITS),(void*)neutral->Z,MAX_BYTES,cudaMemcpyHostToDevice);
		cuda_Memcpy((void*)(iter+3*NB_DIGITS),(void*)neutral->T,MAX_BYTES,cudaMemcpyHostToDevice);
		iter += 4*NB_DIGITS;
	}
	
	// Provést výpočet (sliding window)
	START_MEASURE(start);
	for (int i = coeff.l-1,u,s = 0;i >= 0;)
	{
		if (coeff.bits[i] == 0){
		  curvesDbl<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,swAx);
		  --i;
		}
		else {
			s = i - WINDOW_SZ + 1;
			s = s > 0 ? s : 0;

			while (!coeff.bits[s]) ++s;
			for (int h = 1;h <= i-s+1;++h)  
			  curvesDbl<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,swAx);

			u = coeff.build(s,i);
			if (u > 0){
			  iter = ((digit_t*)swPc)+((u-1)/2)*NUM_CURVES*4*NB_DIGITS;
			  curvesAdd<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,(void*)iter,swAx);
			}
			else { 
			  iter = ((digit_t*)swPc)+((-u-1)/2)*NUM_CURVES*4*NB_DIGITS;
			  curvesSub<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,(void*)iter,swAx); 
			} 
			i = s-1;
		}
	}
	STOP_MEASURE("Computation phase",start,stop,totalTime);
	printf("--------------------------\n");
	printf("Total time: %.3f ms\n",totalTime);

	// Nakopírovat výsledky zpátky do paměti počítače
	iter = (digit_t*)swQw;
	for (int i = 0;i < NUM_CURVES;i++){
	   cuda_Memcpy((void*)initPoints[i].X,(void*)(iter+0*NB_DIGITS),MAX_BYTES,cudaMemcpyDeviceToHost);
	   cuda_Memcpy((void*)initPoints[i].Y,(void*)(iter+1*NB_DIGITS),MAX_BYTES,cudaMemcpyDeviceToHost);
	   cuda_Memcpy((void*)initPoints[i].Z,(void*)(iter+2*NB_DIGITS),MAX_BYTES,cudaMemcpyDeviceToHost);
	   cuda_Memcpy((void*)initPoints[i].T,(void*)(iter+3*NB_DIGITS),MAX_BYTES,cudaMemcpyDeviceToHost);
	   iter += 4*NB_DIGITS;	
	}
 
    // Zkontroluj chyby
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess)
      fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
    
    // Synchronizovat vše
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess)
	  fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
 
	// Uvolnit paměť 
	cuda_Free(swAx);
	cuda_Free(swQw);
	cuda_Free(swPc);

	gpuErrchk(cudaEventDestroy(start));
	gpuErrchk(cudaEventDestroy(stop));

	cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
		return cudaStatus;
    }

    return cudaStatus;
}
