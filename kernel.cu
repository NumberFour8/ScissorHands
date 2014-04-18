#include "kernel.h"

#ifndef USE_TWISTED
	#include "edwards.h"
#else 
	#include "twisted.h"
#endif

// Vrátí výsek z NAF rozvoje
__device__ int build(char* bits,unsigned int start,unsigned int end) 
{
	int ret = 0;
	for (unsigned int i = start;i <= end;i++)
	{
		ret += bits[i]*(1 << (i-start));
	}

	return ret;
}

// Výpočet pomocí sliding window
__global__ void slidingWindow(void* pY,void* pPc,void* swAux,void* swCoeff)
{
	PREPARE();

	VOL digit_t* Qd	  = ((digit_t*)pY)+idx; 
	
	// Nakopírování pracovních dat pro Y
	c_x1[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	char* Cf = (char*)swCoeff;
	for (int i = ax->nafLen-1,u,s = 0;i >= 0;)
	{
		if (Cf[i] == 0)
		{
		  curvesDbl(); 
		  --i;
		}
		else {
			s = i - ax->windowSz + 1;
			s = s > 0 ? s : 0;

			while (!Cf[s]) ++s;
			for (int h = 1;h <= i-s+1;++h) 
			{
			  curvesDbl();
			}

			u = build(Cf,s,i);
			if (u > 0){
			  Qd = ((digit_t*)pPc)+idx+((u-1)/2)*NUM_CURVES*4*NB_DIGITS;
			  curvesAdd();
			}
			else { 
			  Qd = ((digit_t*)pPc)+idx+((-u-1)/2)*NUM_CURVES*4*NB_DIGITS;
			  curvesSub(); 
			} 
			i = s-1;
		}
	}

	// Nakopírování pracovních dat zpátky do Y
	Qd = ((digit_t*)pY) + idx;

	*(Qd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
	*(Qd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
	*(Qd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
	*(Qd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T

	__syncthreads();
}

__global__ void precompute(void* pX,void* pCube,void* swAux)
{
	PREPARE();

	VOL digit_t* Qd    = ((digit_t*)pX)    + idx; 
	VOL digit_t* out   = ((digit_t*)pCube) + idx;

	// Nakopírování pracovních dat pro Y
	c_x1[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	curvesDbl();
	for (int i = 1; i < (1 << (ax->windowSz-1));++i)
	{ 
		curvesAdd();

		// Výsledek na své místo
		*(out+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
		*(out+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
		*(out+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
		*(out+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T
		out += NUM_CURVES*4*NB_DIGITS;
		__syncthreads();

		curvesAdd();
	}
}

cudaError_t compute(const ComputeConfig& h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF& coeff)
{		
	const int PRECOMP_SZ = (1 << (h_input.windowSz-1)); // Počet bodů, které je nutné předpočítat
	const int NUM_CURVES = h_input.numCurves;			// Počet načtených křivek
	const int NUM_BLOCKS = NUM_CURVES/CURVES_BY_BLOCK;	// Počet použitých bloků
	
	cudaEvent_t start,stop;
	float totalTime = 0;
	void *swQw = NULL,*swPc = NULL,*swAx = NULL,*swCf = NULL;
	gpuErrchk(cudaSetDevice(0));
	
	gpuErrchk(cudaEventCreate(&start));
	gpuErrchk(cudaEventCreate(&stop));

	// Alokace potřebných dat
	cuda_Malloc((void**)&swPc,NUM_CURVES*PRECOMP_SZ*4*MAX_BYTES); // Předpočítané body
	cuda_Malloc((void**)&swQw,NUM_CURVES*4*MAX_BYTES);			  // Pomocný bod
	cuda_Malloc((void**)&swAx,sizeof(Aux));						  // Pomocná struktura
	cuda_Malloc((void**)&swCf,h_input.nafLen);					  // NAF rozvoj koeficientu
	
	// Pomocná struktura
	cuda_Memcpy(swAx,(void*)&h_input,sizeof(Aux),cudaMemcpyHostToDevice);

	// NAF rozvoj koeficientu
	cuda_Memcpy(swCf,(void*)coeff.bits,h_input.nafLen,cudaMemcpyHostToDevice);
	
	// Počáteční body
	VOL digit_t* iter = (digit_t*)swPc;
	for (int i = 0;i < NUM_CURVES;i++){
	   cuda_Memcpy((void*)(iter+0*NB_DIGITS),(void*)initPoints[i].X,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+1*NB_DIGITS),(void*)initPoints[i].Y,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+2*NB_DIGITS),(void*)initPoints[i].Z,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+3*NB_DIGITS),(void*)initPoints[i].T,MAX_BYTES,cudaMemcpyHostToDevice);
	   iter += 4*NB_DIGITS;	
	}

	// Konfigurace kernelů
	dim3 threadsPerBlock(NB_DIGITS,CURVES_PER_BLOCK);
	printf("Execution configuration: %d x %d x %d\n",NUM_BLOCKS,CURVES_PER_BLOCK,NB_DIGITS);
	
	// Další předpočítané body
	START_MEASURE(start);
	precompute<<<NUM_BLOCKS,threadsPerBlock>>>((void*)swPc,(void*)iter,(void*)swAx);
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
	slidingWindow<<<NUM_BLOCKS,threadsPerBlock>>>((void*)swQw,(void*)swPc,(void*)swAx,(void*)swCf);
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
	cuda_Free(swCf);

	gpuErrchk(cudaEventDestroy(start));
	gpuErrchk(cudaEventDestroy(stop));

	cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
		return cudaStatus;
    }

    return cudaStatus;
}
