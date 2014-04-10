#include "kernel.h"

// Globální proměnné
__device__ void* d_invN;
__device__ void* d_N;
__device__ void* d_3N;

/*__global__ void edwardsAdd(ExtendedPoint** R, ExtendedPoint **P, ExtendedPoint **Q)
{
	// Proměnné ve sdílené paměti pro bod P
    __shared__ VOL digit_t x1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t y1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t z1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t t1[NUM_CURVES][NB_DIGITS];
	
	// Proměnné ve sdílené paměti pro bod Q
	__shared__ VOL digit_t x2[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t y2[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t z2[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t t2[NUM_CURVES][NB_DIGITS];
	
	// Pomocné proměnné ve sdílené paměti pro přenos a t0,t1,t2
	__shared__ VOL carry_t carry[NUM_CURVES][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t temp1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t temp2[NUM_CURVES][NB_DIGITS];
	
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
	
	VOL digit_t* dd_N = (digit_t*)d_N;
	VOL digit_t* dd_3N = (digit_t*)d_3N;

	const digit_t _N    = dd_N[threadIdx.x];	   // x-tá cifra N
	const digit_t _3N   = dd_3N[threadIdx.x];   // x-tá cifra 3*N
	const digit_t _INVN = *(digit_t*)d_invN;			   // -N^(-1) mod W
	
	// Načítání dat (4 souřadnice po 128 bajtech)
	const digit_t idx = 4*MAX_BYTES*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Teď můžeme přečíst správný bod P
	VOL digit_t* Qd   = ((digit_t*)Q)+idx; // Teď můžeme přečíst správný bod Q
	
	// Nakopírování pracovních dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x);    // prvních 32 cifer patří k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+32); // dalších 32 cifer patří k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+64); // dalších 32 k souřadnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+96); // ... a poslední k souřadnici T
	
	c_tt0[threadIdx.x] = *(Pd+threadIdx.x+32); // t0 = Y1

	c_x2[threadIdx.x] = *(Qd+threadIdx.x);    // prvních 32 cifer patří k X
	c_y2[threadIdx.x] = *(Qd+threadIdx.x+32); // dalších 32 cifer patří k Y
	c_z2[threadIdx.x] = *(Qd+threadIdx.x+64); // dalších 32 k souřadnici Z
	c_t2[threadIdx.x] = *(Qd+threadIdx.x+96); // ... a poslední k souřadnici T
	
	
	c_tcy[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 
	_AUX[threadIdx.x] = 0; 

	// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
	/////////////////////////////////////////	
	
	SUE_MOD(c_tt0,c_x1);
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
	
	// Nakopírování pracovních dat zpátky
	*(Pd+threadIdx.x) 	 = c_x1[threadIdx.x];  // prvních 32 cifer patří k X
	*(Pd+threadIdx.x+32) = c_y1[threadIdx.x];  // dalších 32 cifer patří k Y
	*(Pd+threadIdx.x+64) = c_z1[threadIdx.x];  // dalších 32 k souřadnici Z
	*(Pd+threadIdx.x+96) = c_t1[threadIdx.x];  // ... a poslední k souřadnici T

}*/

__global__ void edwardsDbl(void* R,void* P)
{
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
	
	VOL digit_t* dd_N   = (digit_t*)d_N;
	VOL digit_t* dd_3N  = (digit_t*)d_3N;

	const digit_t _N    = dd_N[threadIdx.x];	// x-tá cifra N
	const digit_t _3N   = dd_3N[threadIdx.x];   // x-tá cifra 3*N
	const digit_t _INVN = *(digit_t*)d_invN;	// -N^(-1) mod W
	
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
/*
void getPrecomputed(** prec,const int exp,ExtendedPoint** pR)
{
	int k = ((exp > 0 ? exp : -exp)-1)/2;
	if (exp > 0)
	 *pR = prec[k];
	//else {
		//mpz_neg(res,precomp[k]);
	//}
}*/

cudaError_t computeExtended(const h_Aux input,h_ExtendedPoint* initPoints,const NAF coeff)
{
	const int PRECOMP_SZ = (1 << (coeff.w-2))+1;		// Počet bodů, které je nutné předpočítat
	const int NUM_CURVES = CURVES_PER_BLOCK*NUM_BLOCKS; // initPoints má tolik prvků
	
	void *swQw = NULL,*swPc = NULL;
	gpuErrchk(cudaSetDevice(0));
	
	cuda_Malloc((void**)swPc,NUM_CURVES*PRECOMP_SZ*4*MAX_BYTES);
	cuda_Malloc((void**)swQw,NUM_CURVES*4*MAX_BYTES);
	
	// Konstanty
	cuda_Malloc((void**)&d_N ,MAX_BYTES);
	cuda_Malloc((void**)&d_3N,MAX_BYTES);
	cuda_Malloc((void**)&d_invN,sizeof(digit_t));
	
	cuda_Memcpy(d_N,(void*)input.N,MAX_BYTES,cudaMemcpyHostToDevice);
	cuda_Memcpy(d_3N,(void*)input.N3,MAX_BYTES,cudaMemcpyHostToDevice);
	cuda_Memcpy(d_invN,(void*)&input.invN,SIZE_DIGIT/8,cudaMemcpyHostToDevice);
	
	
	// Počáteční body
	digit_t* iter = (digit_t*)swPc;
	for (int i = 0;i < NUM_CURVES;i++){
	   cuda_Memcpy((void*)(iter+0*NB_DIGITS),(void*)initPoints->X,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+1*NB_DIGITS),(void*)initPoints->Y,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+2*NB_DIGITS),(void*)initPoints->Z,MAX_BYTES,cudaMemcpyHostToDevice);
	   cuda_Memcpy((void*)(iter+3*NB_DIGITS),(void*)initPoints->T,MAX_BYTES,cudaMemcpyHostToDevice);
	}

	// Další předpočítané body
	iter += 4*NB_DIGITS;
	
	dim3 threadsPerBlock(CURVES_PER_BLOCK, NB_DIGITS);
	edwardsDbl<<<NUM_BLOCKS,threadsPerBlock>>> ((void*)iter,(void*)iter);
	
	/*

    // A počítáme pomocí sliding-window
    int i = coeff.l-1,h,s = 0,k = 0,u;
	while (i >= 0)
	{
		if (coeff.bits[i] == 0){
		  edwardsDbl<<<NUM_CURVES,NB_DIGITS>>>(pts,pts);
		  i--;
		}
		else {
			s = i - coeff.w + 1;
			s = s > 0 ? s : 0;

			while (!coeff.bits[s]) ++s;
			for (h = 1;h <= i-s+1;++h)  
			  edwardsDbl<<<NUM_CURVES,NB_DIGITS>>>(pts,pts);

			u = coeff.build(s,i);

			//getPrecomputed(temp,u);
			//multiply(res,temp);
	
			i = s-1;
		}
	}*/
 
    // Zkontroluj chyby
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess)
      fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
    
    // Synchronizovat vše
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess)
	  fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
 
	cuda_Free(d_N);
	cuda_Free(d_3N);
	cuda_Free(d_invN);

    cuda_Free(initPts);
	cuda_Free(precomp);

    
    // Zkopírovat data zpět do počítače a uvolnit paměť
    for (int i = 0;i < NUM_CURVES;++i){
      pts[i]->toHost(initPoints+i);
      delete pts[i];
    }
    delete[] pts;
    
    return cudaStatus;
}
