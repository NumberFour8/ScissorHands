#ifndef COORDS_H
#define COORDS_H

#include "def.h"
#include "helpers.h"

// Třída pro bod v Extended souřadnicích v paměti počítače
class ExtendedPoint {
private:

	void initAll(bool minusOne);
public:	

	// Je to křivka s a = -1?
	bool isMinus1;

	// Všechny souřadnice bodu v Extended souřadnicích
	biguint_t X,Y,Z,T;
	
	// Vytvoří prázdný bod v Extended souřadnicích
	ExtendedPoint(bool minusOne = true);
	
	// Vytvoří bod v Extended souřadnicích inicializovaný daným afinním bodem
	ExtendedPoint(mpz_t x,mpz_t y,mpz_t N,bool minusOne = true);

	// Vytvoří bod v nekonečnu v Extended souřadnicích
	ExtendedPoint(mpz_t N,bool minusOne = true);
		
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

// Možné strategie výpočtu podle typu načtených křivek
typedef enum { csMixed, csEdwards, csTwisted, csNone } computeStrategy;

/* Nacte krivky a jejich pocatecni body v racionalnich afinnich souradnicich 
 * ze souboru a provede jejich redukci modulo N.
 * Nactene pocatecni body v Extended souradnicich jsou ulozeny do pInit,
 * který je poté nutno uvolnit pomoci delete[].
 * Parametry edwards a twisted jsou počty příslušných typů načtených křivek.
 * Vraci pocet uspesne prectenych krivek. 
*/
computeStrategy readCurves(string file,mpz_t N,ExtendedPoint** pInit,int& edwards,int& twisted,int& usableCurves);


#endif
