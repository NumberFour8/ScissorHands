#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

#include "kernel.h"

using namespace std;

// Přečte počáteční body křivek v afinních souřadnicích ze souboru
int readCurves(string file,mpz_t N,ExtendedPoint** pInit)
{
	ifstream fp;
	if (file.length() < 2) 
	{
		file = "curves_edwards.txt";
		cout << "INFO: Defaulting to curves_edwards.txt" << endl;	
	}
	
	// Pokud se nepodaří otevřít soubor, skonči
	fp.open(file);
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
	mpz_init_set_ui(zX,0);
	mpz_init_set_ui(zY,0);
	
	string ln;
	vector<ExtendedPoint> v;
	bool minus1 = true;
	cout << "Loading curves..." << endl;
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
			if (mpz_cmp_ui(zX,0) != 0) // Byl nalezen faktor?
			{
				cout << "Factor found: " << mpz_to_string(zX) << endl;
			}
			else if (mpz_cmp_ui(zY,0) != 0)
			{
				cout << "Factor found: " << mpz_to_string(zY) << endl;
			}
			
			fp.close();
			mpz_clrs(zX,zY);
			mpq_clrs(qX,qY);

			return 0;
		}

		// Vytvoř bod v Extended souřadnicích z redukovaných afinních bodů modulo N
		v.push_back(ExtendedPoint(zX,zY,N)); 
	}

	// Překroucené Edwardsovy křivky přijdou na začátek
	std::sort(v.begin(), v.end(), [](const ExtendedPoint& a, const ExtendedPoint & b) -> bool { return a.minusOne && !b.minusOne; });

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
	string ln,curveFile;
	typedef pair<string,bool> factor;
	int exitCode = 0,read_curves = 0;
	char c = 0;

	// Množina nalezených faktorů s vlastním uspořádnáním
	set<factor,bool(*)(factor x, factor y)> foundFactors([](factor x, factor y){ return x.first < y.first; });

	// Načíst N
	mpz_t zN;
	cout << "Enter N:" << endl;
	cin >> ln;
	cout << endl;
	mpz_init_set_str(zN,ln.c_str(),10);

	// Inicializace proměnných
	ExtendedPoint infty(zN); // Neutrální prvek
	Aux ax(zN);			     // Pomocná struktura
	NAF S(2);				 // NAF rozvoj šířky 2
	mpz_t zS,zInvW,zX,zY,zF; // Pomocné proměnné
	cudaError_t cudaStatus;	 // Proměnná pro chybové kódy GPU
	ExtendedPoint *PP;		 // Adresa všech bodů

	// Načíst křivky ze souboru
	cout << "Enter path to curve file:" << endl;
	cin >> curveFile;
	cout << endl;

	restart_bound:

	PP = NULL;
	read_curves = readCurves(curveFile,zN,&PP);

	// Zkontroluj počet načtených křivek
	if (read_curves <= 0)
	{
		cout << "ERROR: No curves read." << endl;
		exitCode = 1;
		goto end;
	}
	else if (read_curves != NUM_CURVES)
	{
		cout << "ERROR: Invalid number of curves in file."  
			 << "Required: " << NUM_CURVES << ". Got: " << read_curves 
			 << endl;
		exitCode = 1;
		goto end;
	}
	cout << "Loaded " << read_curves << " curves." << endl << endl;

	// Přečti B1
	cout << "Enter B1: " << endl;
	cin >> ln;
	cout << endl;

	// Spočti S = lcm(1,2,3...,B1) a jeho NAF rozvoj
	mpz_init(zS);
	//mpz_set_ui(zS,(unsigned int)std::stoul(ln));
	
	cout << "Computing coefficient..." << endl;
	lcmToN(zS,(unsigned int)std::stoul(ln));
	S.initialize(zS);
	//cout << "s = " << mpz_to_string(zS) << endl << endl;
	mpz_clear(zS);	
	
	cout << endl << "Computation started..." << endl;
	// Proveď výpočet
	cudaStatus = compute(ax,&infty,PP,S);
    if (cudaStatus != cudaSuccess) 
    {
        fprintf(stderr, "CUDA compute failed!");
		exitCode = 1;
        goto end;
    }
	cout << "Computation finished." << endl;

	// Inicializace pomocných proměnných
	mpz_intz(zInvW,zX,zY,zF);
	
	// Spočti 2^(-W) mod N 
	mpz_ui_pow_ui(zInvW, 2, SIZE_DIGIT); 
	mpz_invert(zInvW, zInvW, zN);

	cout << endl;
	cout << "Print also affine coordinates? (y/n)" << endl;
	cin >> c;
	cout << endl;

	// Analyzuj výsledky
	foundFactors.clear();
	for (int i = 0; i < NUM_CURVES;++i)
	{
		cout << "Curve #" << i+1 << ":\t"; 
		if (PP[i].toAffine(zX,zY,zN,zInvW,zF)) 
		{
			cout << "No factor found." << endl;
			if (c == 'y')
			  cout << endl << "sP = (" << mpz_to_string(zX) << "," << mpz_to_string(zY) << ")" << endl;
		}
		else if (mpz_cmp_ui(zF,0) != 0) // Máme faktor!
		{
		   string fact = mpz_to_string(zF);
		   if (mpz_probab_prime_p(zF,25) != 0)
		   {
		      cout << "Prime factor found: " << fact << endl;
		      foundFactors.insert(std::make_pair(fact,true));
		   }
		   else
		   {
			  cout << "Composite factor found: " << fact << endl;
			  foundFactors.insert(std::make_pair(fact,false));
		   }

		}
		else cout << "Error during conversion." << endl;
		cout << endl << "------------" << endl;
    }
	
	// Vypiš všechny nalezené faktory
	if (foundFactors.size() > 0)
	{
	  cout << endl << foundFactors.size() << " FACTORS FOUND:" << endl << endl;
	  std::for_each(foundFactors.begin(),foundFactors.end(),
				    [](const factor f)
					{ cout << (f.second ? "Prime:\t" : "Composite:\t") << f.first <<  endl; }
				   );
	}
	else cout << endl << "No factors found." << endl << endl;

	// Vyčisti proměnné
	mpz_clrs(zInvW,zX,zY);
	delete[] PP;

	cout << endl << "Type 'r' to restart with new B1 or 'q' to quit..." << endl;
	while (1) {
	   cin >> c;
	   if (c == 'q') break;
	   if (c == 'r') goto restart_bound;
	}
	
	end:
	mpz_clear(zN);

	return exitCode;
}
