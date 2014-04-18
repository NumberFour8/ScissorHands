#ifndef EDWARDS_H
#define EDWARDS_H

#include "def.h"
#include "helpers.h"

// Třída pro bod v Extended souřadnicích v paměti počítače
class ExtendedPoint {
private:

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
	
	// Vytvoří prázdný bod v Extended souřadnicích
	ExtendedPoint()
	{
		initAll();
	}
	
	// Vytvoří bod v Extended souřadnicích inicializovaný daným afinním bodem
	ExtendedPoint(mpz_t x,mpz_t y,mpz_t N) 
	{
		initAll();
		fromAffine(x,y,N);
	}

	// Vytvoří bod v nekonečnu v Extended souřadnicích
	ExtendedPoint(mpz_t N) 
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


	/* Transformace z afinních souřadnic do Extended souřadnic v Montgomeryho reprezentaci
	   Předpoklad: X,Y < N
	*/
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
		V případě chyby vrací false a případný nalezný faktor čísla N je zapsán do fact.
	*/
	bool toAffine(mpz_t x,mpz_t y,mpz_t N,mpz_t invB,mpz_t fact)
	{
		mpz_t z;
		mpz_init(z);
		mpz_set_ui(fact,0);

		biguint_to_mpz(x,X);
		biguint_to_mpz(y,Y);
		biguint_to_mpz(z,Z);

		from_mont_repr(x,N,invB);
		from_mont_repr(y,N,invB);
		from_mont_repr(z,N,invB);
	
		bool ret = try_invert_mod(fact,z,N);
		if (ret) 
		{
			mpz_mul(x,x,fact);
			mpz_mul(y,y,fact);
			mpz_mod(x,x,N);
			mpz_mod(y,y,N);
			
			mpz_set_ui(fact,0);
		}
		
		mpz_clear(z);
		return ret;
	}
};

#endif
