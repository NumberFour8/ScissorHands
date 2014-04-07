#ifndef EDWARDS_H
#define EDWARDS_H

#include "def.h"

// TODO: h_ struktury jako bázové třídy pro GPU třídy

// Struktura pro bod v Inverted souřadnicích v paměti počítače
struct h_InvertedPoint {
	biguint_t X,Y,Z;
	
	void fromMPZ(mpz_t x,mpz_t y,mpz_t z)
	{
		mpz_to_biguint(X,x);
		mpz_to_biguint(Y,y);
		mpz_to_biguint(Z,z);
	}
	
	void toMPZ(mpz_t x,mpz_t y,mpz_t z)
	{
		biguint_to_mpz(x,X);
		biguint_to_mpz(y,Y);
		biguint_to_mpz(z,Z);
	}
};

// Tříd reprezentující bod v Inverted souřadnicích v paměti GPU
class InvertedPoint {
public:
	biguint_t X,Y,Z;
	InvertedPoint() 
	{
		cudaMalloc((void**)&X, MAX_BYTES);
		cudaMalloc((void**)&Y, MAX_BYTES);
		cudaMalloc((void**)&Z, MAX_BYTES);
	}
	
	void toGPU(const h_InvertedPoint* P)
	{
		cudaMemcpy((void*)X,(void*)P->X, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)Y,(void*)P->Y, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)Z,(void*)P->Z, MAX_BYTES, cudaMemcpyHostToDevice);
	}
	
	void toHost(h_InvertedPoint* P) const 
	{
		cudaMemcpy((void*)P->X,(void*)X, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Y,(void*)Y, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Z,(void*)Z, MAX_BYTES, cudaMemcpyDeviceToHost);
	}
	
	virtual ~InvertedPoint()
	{
		cudaFree((void*)X);
		cudaFree((void*)Y);
		cudaFree((void*)Z);
	}
};


// Struktura pro bod v Extended souřadnicích v paměti počítače
struct h_ExtendedPoint {
	biguint_t X,Y,Z,T;
	void fromMPZ(mpz_t x,mpz_t y,mpz_t z,mpz_t t)
	{
		mpz_to_biguint(X,x);
		mpz_to_biguint(Y,y);
		mpz_to_biguint(Z,z);
		mpz_to_biguint(T,t);
	}
	
	void toMPZ(mpz_t x,mpz_t y,mpz_t z,mpz_t t)
	{
		biguint_to_mpz(x,X);
		biguint_to_mpz(y,Y);
		biguint_to_mpz(z,Z);
		biguint_to_mpz(t,T);
	}
};

// Tříd reprezentující bod v Extended souřadnicích v paměti GPU
class ExtendedPoint  {
public:
	biguint_t X,Y,Z,T;
	ExtendedPoint()
	{
		cudaMalloc((void**)&X, MAX_BYTES);
		cudaMalloc((void**)&Y, MAX_BYTES);
		cudaMalloc((void**)&Z, MAX_BYTES);
		cudaMalloc((void**)&T, MAX_BYTES);
	}
	
	void toGPU(const h_ExtendedPoint* P)
	{
		cudaMemcpy((void*)X,(void*)P->X, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)Y,(void*)P->Y, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)Z,(void*)P->Z, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)T,(void*)P->T, MAX_BYTES, cudaMemcpyHostToDevice);
	}
	
	void toHost(h_ExtendedPoint* P) const 
	{
		cudaMemcpy((void*)P->X,(void*)X, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Y,(void*)Y, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Z,(void*)Z, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->T,(void*)T, MAX_BYTES, cudaMemcpyDeviceToHost);
	}
	
	virtual ~ExtendedPoint()
	{
		cudaFree((void*)X);
		cudaFree((void*)Y);
		cudaFree((void*)Z);
		cudaFree((void*)T);
	}
};

#endif
