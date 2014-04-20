#ifndef COORDS_H
#define COORDS_H

#include "def.h"
#include "helpers.h"

// Třída pro bod v Extended souřadnicích v paměti počítače
class ExtendedPoint {
private:

	void initAll();
public:	

	// Všechny souřadnice bodu v Extended souřadnicích
	biguint_t X,Y,Z,T;
	
	// Vytvoří prázdný bod v Extended souřadnicích
	ExtendedPoint();
	
	// Vytvoří bod v Extended souřadnicích inicializovaný daným afinním bodem
	ExtendedPoint(mpz_t x,mpz_t y,mpz_t N);

	// Vytvoří bod v nekonečnu v Extended souřadnicích
	ExtendedPoint(mpz_t N);
		
	// Nastaví na neutrální prvek na Edwardsově křivce
	void infinity(mpz_t N);


	/* Transformace z afinních souřadnic do Extended souřadnic v Montgomeryho reprezentaci
	   Předpoklad: X,Y < N
	*/
	void fromAffine(mpz_t x,mpz_t y,mpz_t N);

	
	/*
		Převede bod z Extended souřadnic v Montgomeryho reprezentaci zpět do afinních.
		V případě chyby vrací false a případný nalezný faktor čísla N je zapsán do fact.
	*/
	bool toAffine(mpz_t x,mpz_t y,mpz_t N,mpz_t invB,mpz_t fact);

};

/* Nacte krivky a jejich pocatecni body v racionalnich afinnich souradnicich 
 * ze souboru a provede jejich redukci modulo N.
 * Nactene pocatecni body v Extended souradnicich jsou ulozeny do pInit,
 * který je poté nutno uvolnit pomoci delete[].
 * Indikator minus1 je nastaven v pripade, ze se jedna o prekroucene Edwardsovy
 * krivky s parametrem a = 1.
 * Vraci pocet uspesne prectenych krivek. 
*/
int readCurves(string file,mpz_t N,ExtendedPoint** pInit,bool& minus1);


#endif
