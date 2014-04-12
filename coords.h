#ifndef EDWARDS_H
#define EDWARDS_H

#include "def.h"
#include "helpers.h"

// Třída pro bod v Extended souřadnicích v paměti počítače
class ExtendedPoint {
private:
	/*	
		Vrací true, je-li v invx platný inverzní prvek k x v okruhu modulo N
		při chybě vrací false a potom invx buď je přímo faktor N, nebo 0.
	*/
	bool inverseMod(mpz_t invx,mpz_t x,mpz_t N)
     {
		 bool ret = false;
		 
		 mpz_t b,r;
		 mpz_init(b);
		 mpz_init(r);
         
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

		 mpz_clear(b);
		 mpz_clear(r);
		 return ret;
     }

	void initAll()
	{
		reset(X);
		reset(Y);
		reset(Z);
		reset(T);
	}
public:	
	// Všechny souřadnice bodu v Extended souřadnicích
	biguint_t X,Y,Z,T;
	
	// Indikátor, zda se počítá na překroucené Edwardsově křivce s a = -1
	bool minusOne;

	// Vytvoří prázdný bod v Extended souřadnicích
	ExtendedPoint(bool minus1 = true)
	  : minusOne(minus1)
	{
		initAll();
	}
	
	// Vytvoří bod v Extended souřadnicích inicializovaný daným afinním bodem
	ExtendedPoint(mpz_t x,mpz_t y,mpz_t N,bool minus1 = true) 
	  : minusOne(minus1)
	{
		initAll();
		fromAffine(x,y,N);
	}
	
	// Vytvoří bod v nekonečnu v Extended souřadnicích
	ExtendedPoint(mpz_t N,bool minus1 = true) 
	  : minusOne(minus1)
	{
		initAll();
		infinity(N);
	}
		
	// Nastaví na neutrální prvek na Edwardsově křivce
	void infinity(mpz_t N)
	{
		mpz_t x,y;
		mpz_init_set_ui(x,0);
		mpz_init_set_ui(y,1);

		fromAffine(x,y,N);
		mpz_clrs(x,y);
	}


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

		mpz_clrs(t,z);
	}

	
	/*
		Převede bod z Extended souřadnic v Montgomeryho reprezentaci zpět do afinních.
		V případě chyby vrací false a případný nalezný faktor N je ve struktuře pRes.
	*/
	bool toAffine(mpz_t x,mpz_t y,mpz_t N,mpz_t invB,ExtResult* pRes)
	{
		mpz_t z,f;
		mpz_intz(z,f);

		pRes->factorFound = false;
		mpz_set_ui(pRes->factor,0);

		biguint_to_mpz(x,X);
		biguint_to_mpz(y,Y);
		biguint_to_mpz(z,Z);

		from_mont_repr(x,N,invB);
		from_mont_repr(y,N,invB);
		from_mont_repr(z,N,invB);
	
		if (!inverseMod(f,z,N)){
			pRes->factorFound = (mpz_cmp_ui(f,0) != 0);
			mpz_set(pRes->factor,f);
			
			mpz_clrs(z,f);
			return false;
		}
		else {
			mpz_mul(x,x,f);
			mpz_mul(y,y,f);
			mpz_mod(x,x,N);
			mpz_mod(y,y,N);

			mpz_clrs(z,f);
			return true;
		}
	}
};

#endif
