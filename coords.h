#ifndef EDWARDS_H
#define EDWARDS_H

#include "def.h"
#include "helpers.h"

// Třída pro bod v Extended souřadnicích v paměti počítače
class h_ExtendedPoint {
private:
	/*	
		Vrací true, je-li v invx platný inverzní prvek k x v okruhu modulo N
		při chybě vrací false a potom invx buď je přímo faktor N, nebo 0.
	*/
	bool inverseMod(mpz_t invx,mpz_t x,mpz_t N)
     {
		 bool ret = false;
		 
		 mpz_t b,r;
		 mpz_inits(b,r);
         
		 mpz_gcdext(r,invx,b,x,N);
         if (mpz_cmp_ui(r,1) == 0){
		   if (mpz_sgn(invx) == -1)
			 mpz_add(invx,invx,N);
		   ret = true;
		 }
         else if (mpz_cmp(r,N) == 0)
		 {
			 mpz_set_si(invx,0);
		 }
		 else mpz_set(invx,r);

		 mpz_clears(b,r);
		 return ret;
     }

public:	
	biguint_t X,Y,Z,T;
	
	// Transformace z afinních souřadnic do Extended souřadnic v Montgomeryho reprezentaci
	void fromAffine(mpz_t x,mpz_t y,mpz_t N)
	{
		mpz_t z,t;
		mpz_init_set_ui(z,1);
		mpz_init(t);
		mpz_mul(t,x,y);
		
		to_mont_repr(x,N);
		to_mont_repr(y,N);
		to_mont_repr(z,N);
		to_mont_repr(t,N);

		mpz_to_biguint(X,x);
		mpz_to_biguint(Y,y);
		mpz_to_biguint(Z,z);
		mpz_to_biguint(T,t);

		mpz_clears(z,t);
	}

	
	/*
		Převede bod z Extended souřadnic v Montgomeryho reprezentaci zpět do afinních.
		V případě chyby vrací false a případný nalezný faktor N je ve struktuře pRes.
	*/
	bool toAffine(mpz_t x,mpz_t y,mpz_t N,mpz_t invB,ExtResult* pRes)
	{
		mpz_t z,f;
		mpz_inits(z,f);

		pRes->factorFound = false;
		mpz_set_ui(pRes->factor,0);

		biguint_to_mpz(x,X);
		biguint_to_mpz(y,Y);
		biguint_to_mpz(z,X);

		from_mont_repr(x,N,invB);
		from_mont_repr(y,N,invB);
		from_mont_repr(z,N,invB);
	
		if (!inverseMod(f,z,N)){
			pRes->factorFound = (mpz_cmp_ui(f,0) != 0);
			mpz_set(pRes->factor,f);
			
			mpz_clears(z,f);
			return false;
		}
		else {
			mpz_mul(x,x,f);
			mpz_mul(y,y,f);
			mpz_mod(x,x,N);
			mpz_mod(y,y,N);

			mpz_clears(z,f);
			return true;
		}

	}
};

// Třída reprezentující bod v Extended souřadnicích v paměti GPU
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
