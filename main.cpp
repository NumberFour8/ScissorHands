#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>

#include "kernel.h"

using namespace std;

// Přečte počáteční body křivek v afinních souřadnicích ze souboru
int readCurves(string file,mpz_t N,ExtendedPoint** pInit)
{
	ifstream fp;
	fp.open(file);

	// Pokud se nepodaří otevřít soubor, skonči
	if (!fp.is_open())
	{
		cout << "ERROR: Cannot open file." << endl;
		return 0;
	} 
	
	// Inicializace racionálních souřadnic
	mpq_t qX,qY;
	mpq_intz(qX,qY);

	// Inicializace celočíselných souřadnic
	mpz_t zX,zY;
	mpz_intz(zX,zY);
	
	string ln;
	vector<ExtendedPoint> v;
	bool minus1 = true;
	while (getline(fp,ln))
	{
		// Přeskoč segment, který nezačíná #
		if (ln.find("#") == string::npos) continue;
		
		// Je to překroucená Edwardsova křivka s a = -1 ? 
		fp >> ln; 
		minus1 = (ln == "-1");
		
		// Přečti racionální X-ovou souřadnici a zkrať
		fp >> ln;
		mpq_set_str(qX,ln.c_str(),10);
		mpq_canonicalize(qX);
		
		// Přečti racionální Y-ovou souřadnici a zkrať
		fp >> ln;
		mpq_set_str(qY,ln.c_str(),10);
		mpq_canonicalize(qY);

		// Pokus se X-ovou a Y-ovou souřadnici rekudovat modulo N
		if (!reduce_mod(zX,qX,N) || !reduce_mod(zY,qY,N))
		{
			cout << "ERROR: Cannot reduce on curve #" << v.size() << endl;
			fp.close();
			mpz_clrs(zX,zY);
			mpq_clrs(qX,qY);

			return 0;
		}

		// Vytvoř bod v Extended souřadnicích z redukovaných afinních bodů modulo N
		v.push_back(ExtendedPoint(zX,zY,N)); 
	}

	// Překopíruj body z vektoru do paměti
	*pInit = new ExtendedPoint[v.size()];
	std::copy(v.begin(),v.end(),*pInit);

	// Vyčisti paměť
	fp.close();
	mpz_clrs(zX,zY);
	mpq_clrs(qX,qY);

	return (int)v.size();
}

int main()
{
	string ln;

	// Načíst N
	mpz_t zN;
	cout << "Enter N:" << endl;
	cin >> ln;
	mpz_init_set_str(zN,ln.c_str(),10);

	
	// Načíst křivky ze souboru
	cout << "Enter path to curve file:" << endl;
	cin >> ln;

	ExtendedPoint *PP = NULL;
	int read_curves = readCurves(ln,zN,&PP);

	// Zkontroluj počet načtených křivek
	if (read_curves <= 0)
	{
		cout << "ERROR: No curves read." << endl;
		return 1;
	}
	else if (read_curves != NUM_CURVES)
	{
		cout << "ERROR: Invalid number of curves in file."  
			 << "Required: " << NUM_CURVES << ". Got: " << read_curves 
			 << endl;
		return 1;
	}
	
	// Inicializace B1
	mpz_t zS;
	mpz_init(zS);
	
	// Přečti B1 a spočti S = lcm(1,2,3...,B1)
	cout << "Enter B1:" << endl;
	cin >> ln;
	lcmToN(zS,(unsigned int)std::stoul(ln));
	
	// ... a vypočítej NAF rozvoj čísla S
	NAF S(2,zS);
	cout << "S = ";
	printmpz("%s\n",zS);
	mpz_clear(zS);	
		
	// Neturální prvek a pomocná struktura	
	ExtendedPoint infty(zN);
	Aux ax(zN);
	
	// Proveď výpočet
	cudaError_t cudaStatus = compute(ax,&infty,PP,S);
    if (cudaStatus != cudaSuccess) 
    {
        fprintf(stderr, "CUDA compute failed!");
        return 1;
    }

	// Inicializace pomocných proměnných
	mpz_t zInvW,zX,zY;
	mpz_intz(zInvW,zX,zY);
	
	// Spočti 2^(-W) mod N 
	mpz_ui_pow_ui(zInvW, 2, SIZE_DIGIT); 
	mpz_invert(zInvW, zInvW, zN);

	cout << endl;
	
	// Analyzuj výsledky
	for (int i = 0; i < NUM_CURVES;++i)
	{
		cout << "Result at curve #" << i+1 << ":" << endl;
		if (PP[i].toAffine(zX,zY,zN,zInvW))
		{
			printmpz("Affine point: (%s,",zX);
			printmpz("%s)\n",zY);
		}
		else break; // Převod do afinních souřadnic selhal, máme faktor
		cout << endl;
    }
	
	// Vyčisti paměť
	mpz_clrs(zN,zInvW,zX,zY);
	delete[] PP;

	cout << "Press Enter to quit..." << endl;
	cin.ignore();
	
    return 0;
}
