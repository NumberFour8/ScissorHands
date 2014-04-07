#include "kernel.h"

// Globální proměnné
__constant__ __device__ digit_t d_invN;
__device__ biguint_t d_N;
__device__ biguint_t d_3N;


__global__ void edwardsAdd(ExtendedPoint* R, ExtendedPoint *P, ExtendedPoint *Q)
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
	
	VOL digit_t* _CARRY = carry[threadIdx.y];  // přenos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocná proměnná pro násobení
	
	const digit_t _N    = d_N[threadIdx.x];	   // x-tá cifra N
	const digit_t _3N   = d_3N[threadIdx.x];   // x-tá cifra 3*N
	const digit_t _INVN = d_invN;			   // -N^(-1) mod W
	
	// Nakopírování pracovních dat
	const digit_t idx = blockIdx.x*blockDim.y + threadIdx.y;
	
	c_x1[threadIdx.x] = P[idx].C.X[threadIdx.x];
	c_y1[threadIdx.x] = P[idx].C.Y[threadIdx.x];
	c_z1[threadIdx.x] = P[idx].C.Z[threadIdx.x];
	c_t1[threadIdx.x] = P[idx].C.T[threadIdx.x];

	c_x2[threadIdx.x] = Q[idx].C.X[threadIdx.x];
	c_y2[threadIdx.x] = Q[idx].C.Y[threadIdx.x];
	c_z2[threadIdx.x] = Q[idx].C.Z[threadIdx.x];
	c_t2[threadIdx.x] = Q[idx].C.T[threadIdx.x];

	c_cy[threadIdx.x] = 0;
	c_t0[threadIdx.x] = P[idx].y[threadIdx.x]; // t0 = Y1
	c_t1[threadIdx.x] = 0; 
	c_t2[threadIdx.x] = 0; 

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
	ADD_MOD(c_t2,c_tt1,c_tt2);
	
	MUL_MOD(c_x1,c_y2,c_z2);
	MUL_MOD(c_y1,c_t2,c_x2);
	MUL_MOD(c_t1,c_y2,c_x2);
	MUL_MOD(c_z1,c_z2,c_t2);
	
	/////////////////////////////////////////
	R[idx].C.X[threadIdx.x] = c_x1[threadIdx.x];
	R[idx].C.Y[threadIdx.x] = c_y1[threadIdx.x];
	R[idx].C.Z[threadIdx.x] = c_z1[threadIdx.x];
	R[idx].C.T[threadIdx.x] = c_t1[threadIdx.x];
}

__global__ void edwardsDbl(ExtendedPoint*R, ExtendedPoint *P)
{
    // Proměnné ve sdílené paměti pro bod P
    __shared__ VOL digit_t x1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t y1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t z1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t t1[NUM_CURVES][NB_DIGITS];
	
	// Pomocné proměnné ve sdílené paměti pro přenos a t0,t1,t2
	__shared__ VOL carry_t carry[NUM_CURVES][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t temp1[NUM_CURVES][NB_DIGITS];
	__shared__ VOL digit_t temp2[NUM_CURVES][NB_DIGITS];
	
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];
		
	// Pomocné proměnné a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	
	VOL digit_t* _CARRY = carry[threadIdx.y];  // přenos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocná proměnná pro násobení
	
	const digit_t _N    = d_N[threadIdx.x];	   // x-tá cifra N
	const digit_t _3N   = d_3N[threadIdx.x];   // x-tá cifra 3*N
	const digit_t _INVN = d_invN;			   // -N^(-1) mod W
	
	// Nakopírování pracovních dat	
	c_x1[threadIdx.x] = P[idx].C.X[threadIdx.x];
	c_y1[threadIdx.x] = P[idx].C.Y[threadIdx.x];
	c_z1[threadIdx.x] = P[idx].C.Z[threadIdx.x];
	c_t1[threadIdx.x] = P[idx].C.T[threadIdx.x];

	c_cy[threadIdx.x] = 0;
	c_t0[threadIdx.x] = 0;
	c_t1[threadIdx.x] = 0; 
	c_t2[threadIdx.x] = 0; 
 
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
	R[idx].C.X[threadIdx.x] = c_x1[threadIdx.x];
	R[idx].C.Y[threadIdx.x] = c_y1[threadIdx.x];
	R[idx].C.Z[threadIdx.x] = c_z1[threadIdx.x];
	R[idx].C.T[threadIdx.x] = c_t1[threadIdx.x];
}

void aux_getPointMultiples(ExtendedPoint* R,ExtendedPoint *P,const unsigned int multiple)
{
	edwardsDbl<<NUM_CURVES,NB_DIGITS>>(R,P);
	if (multiple == 2) return;
	for (int i = 3;i <= multiple;++i){
	  edwardsAdd<<NUM_CURVES,NB_DIGITS>>(R,R,P);
	}
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

void getPrecomputed(const ExtendedPoint** prec,const int exp,ExtendedPoint** pR)
{
	int k = ((exp > 0 ? exp : -exp)-1)/2;
	if (exp > 0){
	 *pR = prec[k];
	//else {
		//mpz_neg(res,precomp[k]);
	//}
}

extern "C" cudaError_t computeExtended(const h_Aux input,h_ExtendedPoint* initPoints,const NAF coeff)
{
	cudaError_t cudaStatus;
    
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return cudaStatus;
    }

	// Konstanty do konstantní paměti
	cudaMemcpyToSymbol(d_N,   (void*)input.N, MAX_BYTES);
	cudaMemcpyToSymbol(d_3N,  (void*)input.N3,MAX_BYTES);
	cudaMemcpyToSymbol(d_invN,(void*)input.invN, SIZE_DIGIT/8);

	// Nakopírovat výchozí body do paměti GPU
	ExtendedPoint **pts = new ExtendedPoint[NUM_CURVES];
	for (int i = 0;i < NUM_CURVES;++i){
	   pts[i] = new ExtendedPoint();
	   pts[i].toGPU(initPoints[i]);
	}
    
    // Předpočítat body pro sliding window
    int precompSize = (1 << (coeff.w-2))+1;
    ExtendedPoint **prec = new ExtendedPoint[precompSize*NUM_CURVES];
    for (int i = 0; i < precompSize;++i){
	   prec[i] = new ExtendedPoint();
	   aux_getPointMultiples(prec[i],pts,2*i+1);
	}
    
    // A počítáme pomocí sliding-window
    int i = exp.length-1,h,s = 0,k = 0,u;
	while (i >= 0)
	{
		if (exp.bits[i] == 0){
		  edwardsDbl(P,P);
		  i--;
		}
		else {
			s = i - w + 1;
			s = s > 0 ? s : 0;

			while (!exp.bits[s]) ++s;
			for (h = 1;h <= i-s+1;++h) square(res);

			u = buildFromNAF(coeff,s,i);

			getPrecomputed(temp,u);
			multiply(res,temp);
			counter++;
			i = s-1;
		}
	}
    
        
    
    // Zkontroluj chyby
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess)
      fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
    
    // Synchronizovat vše
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess)
	  fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
 
    // Zkopírovat data zpět do počítače a uvolnit paměť
    for (int i = 0;i < NUM_CURVES;++i){
      pts[i].toHost(initPoints+i);
      delete pts[i];
    }
    delete[] pts;
    
    return cudaStatus;
}
