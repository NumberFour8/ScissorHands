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
	fp.open(ln);

	if (!fp.is_open())
	{
		cout << "ERROR: Cannot open file." << endl;
		return 0;
	} 
	
	mpz_t zX,zY;
	mpz_intz(zX,zY);
	
	vector<ExtendedPoint> v;
	while (getline(ln,fp))
	{
		if (ln.find("#") == string::npos) continue;
		
		fp >> X;
		fp >> Y;
		mpz_set_str(zX,X.c_str(),10);
		mpz_set_str(zY,Y.c_str(),10);
		
		v.push_back(ExtendedPoint(zX,zY,N)); 
	}
	fp.close();
	mpz_clrs(zX,zY);

	*pInit = new ExtendedPoint[v.size()];
	std::copy(v.begin(),v.end(),*pInit);
	
	//std::for_each(v.begin(),v.end(), [](ExtendedPoint p) { delete &p });
	
	return (int)v.size();
}

int main()
{
	string ln;
	cout << "Enter N:" << endl;
	cin >> ln;

	mpz_t zN;
	mpz_init_set_str(zN,ln.c_str(),10);

	cout << "Enter path to curve file:" << endl;
	cin >> ln;

	ExtendedPoint *PP = NULL;
	
	int read_curves = readCurves(ln,zN,&PP);
	if (read_curves <= 0)
	{
		cout << "ERROR: No curves read" << endl;
		return 1;
	}
	else if (read_curves != NUM_CURVES)
	{
		cout << "ERROR: Invalid number of curves in file."  
			 << "Required: " << NUM_CURVES << ". Got: " << read_curves 
			 << endl;
		return 1;
	}
	
	mpz_t zS;
	mpz_init(zS);
	
	// Hranice první fáze
	cout << "Enter B1:" << endl;
	cin >> ln;
	lcmToN(zS,(unsigned int)std::stoul(ln));
	
	// ... a NAF rozvoj S = lcm(1,2,3...B1)
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

	mpz_t zInvW,zX,zY;
	mpz_intz(zInvW,zX,zY);
	
	// 2^(-W) mod N 
	mpz_ui_pow_ui(zW, 2, SIZE_DIGIT); 
	mpz_invert(zInvW, zInvW, zN);

	cout << endl;
	
	// Analyzuj výsledky
	ExtResult Result;
	for (int i = 0; i < NUM_CURVES;++i)
	{
		if (PP[i].toAffine(zX,zY,zN,zInvW,&Result)){
			cout << "No factor found by curve #" << i+1 << endl;
			printmpz("Affine point: (%s,",zX);
			printmpz("%s)\n",zY);
		}
		else if (Result.factorFound) 
		{
			printmpz("Found factor %s ",res.factor);
			cout << "using curve #" << i+1 << endl;
			break;
		}
		else cout << "ERROR while working with curve #" << i+1 << endl;
    }
	
	// Vyčisti paměť
	mpz_clrs(zN,zInvW,zX,zY);
	delete[] PP;

	cout << "Press Enter to quit..." << endl;
	cin.ignore();
	
    return 0;
}
