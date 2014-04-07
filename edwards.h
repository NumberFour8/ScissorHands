#ifndef EDWARDS_H
#define EDWARDS_H

#include "def.h"

// TODO: h_ struktury jako bázové třídy pro GPU třídy

// Struktura pro bod v Inverted souřadnicích v paměti počítače
struct h_InvertedPoint {
	biguint_t X,Y,Z;
};

// Tříd reprezentující bod v Inverted souřadnicích v paměti GPU
class InvertedPoint {
public:
	h_InvertedPoint C;
	InvertedPoint() 
	{
		cudaMalloc((void**)&C.X, MAX_BYTES);
		cudaMalloc((void**)&C.Y, MAX_BYTES);
		cudaMalloc((void**)&C.Z, MAX_BYTES);
	}
	
	void toGPU(const h_InvertedPoint* P)
	{
		cudaMemcpy((void*)C.X,(void*)P->X, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)C.Y,(void*)P->Y, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)C.Z,(void*)P->Z, MAX_BYTES, cudaMemcpyHostToDevice);
	}
	
	void toHost(h_InvertedPoint* P) const 
	{
		cudaMemcpy((void*)P->X,(void*)C.X, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Y,(void*)C.Y, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Z,(void*)C.Z, MAX_BYTES, cudaMemcpyDeviceToHost);
	}
	
	virtual ~InvertedPoint()
	{
		cudaFree((void*)C.X);
		cudaFree((void*)C.Y);
		cudaFree((void*)C.Z);
	}
};


// Struktura pro bod v Extended souřadnicích v paměti počítače
struct h_ExtendedPoint {
	bigint_t X,Y,Z,T;
};

// Tříd reprezentující bod v Inverted souřadnicích v paměti GPU
class ExtendedPoint  {
public:
	h_ExtendedPoint C;
	ExtendedPoint()
	{
		cudaMalloc((void**)&C.X, MAX_BYTES);
		cudaMalloc((void**)&C.Y, MAX_BYTES);
		cudaMalloc((void**)&C.Z, MAX_BYTES);
		cudaMalloc((void**)&C.T, MAX_BYTES);
	}
	
	void toGPU(const h_ExtendedPoint* P)
	{
		cudaMemcpy((void*)C.X,(void*)P->X, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)C.Y,(void*)P->Y, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)C.Z,(void*)P->Z, MAX_BYTES, cudaMemcpyHostToDevice);
		cudaMemcpy((void*)C.T,(void*)P->T, MAX_BYTES, cudaMemcpyHostToDevice);
	}
	
	void toHost(h_ExtendedPoint* P) const 
	{
		cudaMemcpy((void*)P->X,(void*)C.X, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Y,(void*)C.Y, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->Z,(void*)C.Z, MAX_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)P->T,(void*)C.T, MAX_BYTES, cudaMemcpyDeviceToHost);
	}
	
	virtual ~ExtendedPoint()
	{
		cudaFree((void*)C.X);
		cudaFree((void*)C.Y);
		cudaFree((void*)C.Z);
		cudaFree((void*)C.T);
	}
};

#endif
