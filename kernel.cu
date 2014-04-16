#include "kernel.h"

#ifndef USE_TWISTED
	#include "edwards.h"
#else 
	#include "twisted.h"
#endif

#define PREPARE() \
	Aux *ax = (Aux*)swAux;												\
	__shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t x2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t y2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t z2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL digit_t t2[CURVES_PER_BLOCK][NB_DIGITS];				\
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS];			\
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];			\
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];			\
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];			\
	VOL digit_t* c_x1 = x1[threadIdx.y];								\
	VOL digit_t* c_y1 = y1[threadIdx.y];								\
	VOL digit_t* c_z1 = z1[threadIdx.y];								\
	VOL digit_t* c_t1 = t1[threadIdx.y];								\
	VOL digit_t* c_x2 = x2[threadIdx.y];								\
	VOL digit_t* c_y2 = y2[threadIdx.y];								\
	VOL digit_t* c_z2 = z2[threadIdx.y];								\
	VOL digit_t* c_t2 = t2[threadIdx.y];								\
	VOL digit_t* c_tt0  = temp0[threadIdx.y];							\
	VOL digit_t* c_tt1  = temp1[threadIdx.y];							\
	VOL carry_t* _CARRY = carry[threadIdx.y];							\
	VOL digit_t* _AUX   = temp2[threadIdx.y];							\
	const digit_t _N    = ax->N[threadIdx.x];							\
	const digit_t _3N   = ax->N3[threadIdx.x];							\
	const digit_t _INVN = ax->invN;										\
	const digit_t idx = 4*NB_DIGITS*(blockIdx.x*blockDim.y + threadIdx.y);


__global__ void slidingWindow(void* pY,void* pPc,void* swAux,const int WS,const NAF& coeff)
{
	PREPARE();

	VOL digit_t* Qd	  = ((digit_t*)pY)+idx; 

	// Nakopirovani pracovnich dat pro Y
	c_x1[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	for (int i = coeff.l-1,u,s = 0;i >= 0;)
	{
		if (coeff.bits[i] == 0)
		{
		  curvesDbl(); 
		  --i;
		}
		else {
			s = i - WS + 1;
			s = s > 0 ? s : 0;

			while (!coeff.bits[s]) ++s;
			for (int h = 1;h <= i-s+1;++h) 
			{
			  curvesDbl();
			}

			u = coeff.build(s,i);
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

__global__ void precompute(void* pX,void* pCube,void* swAux,const int PS)
{
	PREPARE();

	VOL digit_t* Qd    = ((digit_t*)pX)    + idx; 
	VOL digit_t* out   = ((digit_t*)pCube) + idx;

	// Nakopirovani pracovnich dat pro Y
	c_x1[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	curvesDbl();
	for (int i = 1; i < PS;++i)
	{ 
		curvesAdd();

		// Vysledek na sve misto
		*(out+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
		*(out+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
		*(out+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
		*(out+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T
		out += NUM_CURVES*4*NB_DIGITS;
		__syncthreads();

		curvesAdd();
	}
}

cudaError_t compute(const Aux h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF& coeff,const unsigned int WS)
{
	const int PRECOMP_SZ = (1 << WS)-1;			// Počet bodů, které je nutné předpočítat
	
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
	precompute<<<NUM_BLOCKS,threadsPerBlock>>>((void*)swPc,(void*)iter,(void*)swAx,PRECOMP_SZ);
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
	slidingWindow<<<NUM_BLOCKS,threadsPerBlock>>>((void*)swQw,(void*)swPc,(void*)swAx,WS,&coeff);
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
