#include "kernel.h"

// Globální proměnné
__device__ __constant__ unsigned int d_invN;
__device__ unsigned int* d_N;
__device__ unsigned int* d_3N;

__global__ void edwardsAdd(void* R, void *P, void *Q,void* aux)
{
	Aux *ax = (Aux*)aux;
	
	// Proměnné ve sdílené paměti pro bod P
    __shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];

	__shared__ VOL digit_t x2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Pomocné proměnné ve sdílené paměti pro přenos a t0,t1,t2
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Výsek pro konkrétní křivku
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];

	VOL digit_t* c_x2 = x2[threadIdx.y];
	VOL digit_t* c_y2 = y2[threadIdx.y];
	VOL digit_t* c_z2 = z2[threadIdx.y];
	VOL digit_t* c_t2 = t2[threadIdx.y];
		
	// Pomocné proměnné a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	VOL carry_t* _CARRY = carry[threadIdx.y];  // přenos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocná proměnná pro násobení
	
	const digit_t _N    = ax->N[threadIdx.x];	// x-tá cifra N
	const digit_t _3N   = ax->N3[threadIdx.x];  // x-tá cifra 3*N
	const digit_t _INVN = ax->invN;				// -N^(-1) mod W
	
	// Načítání dat (4 souřadnice po MAX_BYTES bajtech)
	const digit_t idx = 4*MAX_BYTES*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Teď můžeme přečíst správný bod P
	VOL digit_t* Qd   = ((digit_t*)Q)+idx; // Teď můžeme přečíst správný bod P


	// Nakopírování pracovních dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	c_x2[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y2[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z2[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t2[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	c_tcy[threadIdx.x] = 0;
	c_tt0[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 
	_AUX[threadIdx.x]  = 0; 

	// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
	/////////////////////////////////////////	
	
	SUB_MOD(c_tt0,c_y1,c_x1);
	ADD_MOD(c_tt1,c_y2,c_x2);
	
	MUL_MOD(c_tt0,c_tt0,c_tt1);
	ADD_MOD(c_tt1,c_y1,c_x1);
	
	SUB_MOD(c_x1,c_y2,c_x2);
	MUL_MOD(c_tt1,c_tt1,c_x1);
	
	DBL_MOD(c_z2);
	DBL_MOD(c_t2);
	
	MUL_MOD(c_z1,c_z1,c_t2);
	MUL_MOD(c_z2,c_z2,c_t1);
	
	ADD_MOD(c_y2,c_z2,c_z1);
	SUB_MOD(c_x2,c_z2,c_z1);
	
	SUB_MOD(c_z2,c_tt1,c_tt0);
	ADD_MOD(c_t2,c_tt1,c_tt0);
	
	MUL_MOD(c_x1,c_y2,c_z2);
	MUL_MOD(c_y1,c_t2,c_x2);
	MUL_MOD(c_t1,c_y2,c_x2);
	MUL_MOD(c_z1,c_z2,c_t2);
	
	/////////////////////////////////////////
	
	VOL digit_t* Rd   = ((digit_t*)R)+idx; // Teď můžeme přečíst správný bod R

	// Nakopírování pracovních dat zpátky
	*(Rd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
	*(Rd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
	*(Rd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
	*(Rd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T

	__syncthreads();
}

__global__ void edwardsSub(void* R, void *P, void *Q,void* aux)
{
	Aux *ax = (Aux*)aux;
	
	// Proměnné ve sdílené paměti pro bod P
    __shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];

	__shared__ VOL digit_t x2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Pomocné proměnné ve sdílené paměti pro přenos a t0,t1,t2
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Výsek pro konkrétní křivku
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];

	VOL digit_t* c_x2 = x2[threadIdx.y];
	VOL digit_t* c_y2 = y2[threadIdx.y];
	VOL digit_t* c_z2 = z2[threadIdx.y];
	VOL digit_t* c_t2 = t2[threadIdx.y];
		
	// Pomocné proměnné a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	VOL carry_t* _CARRY = carry[threadIdx.y];  // přenos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocná proměnná pro násobení
	
	const digit_t _N    = ax->N[threadIdx.x];	// x-tá cifra N
	const digit_t _3N   = ax->N3[threadIdx.x];  // x-tá cifra 3*N
	const digit_t _INVN = ax->invN;				// -N^(-1) mod W
	
	// Načítání dat (4 souřadnice po MAX_BYTES bajtech)
	const digit_t idx = 4*MAX_BYTES*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Teď můžeme přečíst správný bod P
	VOL digit_t* Qd   = ((digit_t*)Q)+idx; // Teď můžeme přečíst správný bod P


	// Nakopírování pracovních dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	c_x2[threadIdx.x] = ax->N[threadIdx.x];
	c_t2[threadIdx.x] = ax->N[threadIdx.x];

	c_y2[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z2[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z

	c_tcy[threadIdx.x] = 0;
	c_tt0[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS);;
	c_tt1[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); 
	_AUX[threadIdx.x]  = 0; 

	// Potřebujeme vyrobit -X a -T pro bod Q
	SUE_MOD(c_x2,c_tt0); // X = N-X
	SUE_MOD(c_t2,c_tt1); // T = N-T

	c_tt0[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 

	// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
	/////////////////////////////////////////	
	
	SUB_MOD(c_tt0,c_y1,c_x1);
	ADD_MOD(c_tt1,c_y2,c_x2);
	
	MUL_MOD(c_tt0,c_tt0,c_tt1);
	ADD_MOD(c_tt1,c_y1,c_x1);
	
	SUB_MOD(c_x1,c_y2,c_x2);
	MUL_MOD(c_tt1,c_tt1,c_x1);
	
	DBL_MOD(c_z2);
	DBL_MOD(c_t2);
	
	MUL_MOD(c_z1,c_z1,c_t2);
	MUL_MOD(c_z2,c_z2,c_t1);
	
	ADD_MOD(c_y2,c_z2,c_z1);
	SUB_MOD(c_x2,c_z2,c_z1);
	
	SUB_MOD(c_z2,c_tt1,c_tt0);
	ADD_MOD(c_t2,c_tt1,c_tt0);
	
	MUL_MOD(c_x1,c_y2,c_z2);
	MUL_MOD(c_y1,c_t2,c_x2);
	MUL_MOD(c_t1,c_y2,c_x2);
	MUL_MOD(c_z1,c_z2,c_t2);
	
	/////////////////////////////////////////
	
	VOL digit_t* Rd   = ((digit_t*)R)+idx; // Teď můžeme přečíst správný bod R

	// Nakopírování pracovních dat zpátky
	*(Rd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
	*(Rd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
	*(Rd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
	*(Rd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T

	__syncthreads();
}

__global__ void edwardsDbl(void* R,void* P,void* aux)
{
	Aux *ax = (Aux*)aux;
	
	// Proměnné ve sdílené paměti pro bod P
    __shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Pomocné proměnné ve sdílené paměti pro přenos a t0,t1,t2
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Výsek pro konkrétní křivku
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];
		
	// Pomocné proměnné a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	VOL carry_t* _CARRY = carry[threadIdx.y];  // přenos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocná proměnná pro násobení
	
	const digit_t _N    = ax->N[threadIdx.x];	// x-tá cifra N
	const digit_t _3N   = ax->N3[threadIdx.x];  // x-tá cifra 3*N
	const digit_t _INVN = ax->invN;				// -N^(-1) mod W
	
	// Načítání dat (4 souřadnice po MAX_BYTES bajtech)
	const digit_t idx = 4*MAX_BYTES*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Teď můžeme přečíst správný bod P

	
	// Nakopírování pracovních dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x+0*NB_DIGITS); // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+1*NB_DIGITS); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+2*NB_DIGITS); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+3*NB_DIGITS); // ... a poslední k souřadnici T

	c_tcy[threadIdx.x] = 0;
	c_tt0[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 
	_AUX[threadIdx.x]  = 0; 
 
	// Twisted Edwards Extended (dbl-2008-hwcd-4), a = -1, independent of d,incomplete
	/////////////////////////////////////////
	
	ADD_MOD(c_tt0,c_x1,c_y1);
	SQR_MOD(c_tt1,c_x1);

	SQR_MOD(c_x1,c_y1);
	SQR_MOD(c_y1,c_z1);

	ADD_MOD(c_t1,c_tt1,c_x1);
	SUB_MOD(c_z1,c_tt1,c_x1);

	SQR_MOD(c_tt1,c_tt0);
	DBL_MOD(c_y1);

	SUB_MOD(c_tt0,c_t1,c_tt1);
	ADD_MOD(c_tt1,c_y1,c_z1);

	MUL_MOD(c_x1,c_tt1,c_tt0);
	MUL_MOD(c_y1,c_t1,c_z1);
	MUL_MOD(c_t1,c_t1,c_tt0);
	MUL_MOD(c_z1,c_z1,c_tt1);
	
	////////////////////////////////////////

	VOL digit_t* Rd   = ((digit_t*)R)+idx; // Teď můžeme přečíst správný bod R

	// Nakopírování pracovních dat zpátky
	*(Rd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
	*(Rd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
	*(Rd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
	*(Rd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T

	__syncthreads();
}


int buildFromNAF(NAF N,int start,int end)
{
	int i,ret = 0;
	for (i = start;i <= end;i++)
	{
		ret += N.bits[i]*(1 << (i-start));
	}

	return ret;
}


cudaError_t compute(const Aux h_input,const ExtendedPoint* neutral,ExtendedPoint* initPoints,const NAF coeff)
{
	const int WINDOW_SZ	 = 4;							// Velikost okna
	const int PRECOMP_SZ = (1 << (WINDOW_SZ-2))+1;		// Počet bodů, které je nutné předpočítat
	const int NUM_CURVES = CURVES_PER_BLOCK*NUM_BLOCKS; // initPoints má tolik prvků
	
	void *swQw = NULL,*swPc = NULL,*swAx = NULL;
	gpuErrchk(cudaSetDevice(0));
	
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

	printBigInt("X",initPoints[0].X);
	printBigInt("Y",initPoints[0].Y);
	printBigInt("Z",initPoints[0].Z);
	printBigInt("T",initPoints[0].T);

	// Další předpočítané body
	dim3 threadsPerBlock(NB_DIGITS,CURVES_PER_BLOCK);
	edwardsDbl<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)swQw,(void*)swPc,(void*)swAx);
	for (int i = 1; i < PRECOMP_SZ;++i){ // Tady už je iter nastavené na pozici prvních lichých mocnin
		edwardsAdd<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)iter,(void*)swQw,(void*)swPc,(void*)swAx); 
		edwardsAdd<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)swQw,(void*)iter,(void*)swPc,(void*)swAx);
		iter += NUM_CURVES*4*NB_DIGITS;
	} 
	
	// Do swQw nakopírovat neutrální prvek
	iter = (digit_t*)swQw;
	cuda_Memcpy((void*)(iter+0*NB_DIGITS),(void*)neutral->X,MAX_BYTES,cudaMemcpyHostToDevice);
	cuda_Memcpy((void*)(iter+1*NB_DIGITS),(void*)neutral->Y,MAX_BYTES,cudaMemcpyHostToDevice);
	cuda_Memcpy((void*)(iter+2*NB_DIGITS),(void*)neutral->Z,MAX_BYTES,cudaMemcpyHostToDevice);
	cuda_Memcpy((void*)(iter+3*NB_DIGITS),(void*)neutral->T,MAX_BYTES,cudaMemcpyHostToDevice);

	// Provést výpočet (sliding window)
	for (int i = coeff.l-1,u,s;i >= 0;i = s-1)
	{
		if (coeff.bits[i] == 0){
		  edwardsDbl<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,swAx);
		  --i;
		}
		else {
			s = i - WINDOW_SZ + 1;
			s = s > 0 ? s : 0;

			while (!coeff.bits[s]) ++s;
			for (int h = 1;h <= i-s+1;++h)  
			  edwardsDbl<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,swAx);

			u = coeff.build(s,i);
			if (u > 0){
			  iter = (digit_t*)swPc+((u-1)/2);
			  edwardsAdd<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,(void*)iter,swAx);
			}
			else { 
			  iter = (digit_t*)swPc+((-u-1)/2);
			  edwardsSub<<<NUM_BLOCKS,threadsPerBlock>>>(swQw,swQw,(void*)iter,swAx); 
			} 
		}
	}

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

	cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
		return cudaStatus;
    }

    return cudaStatus;
}
