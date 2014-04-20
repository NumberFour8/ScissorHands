#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#include "coords.h"

void ExtendedPoint::initAll()
{
	reset(X);
	reset(Y);
	reset(Z);
	reset(T);
}

ExtendedPoint::ExtendedPoint()
{
	initAll();
}
	
ExtendedPoint::ExtendedPoint(mpz_t x,mpz_t y,mpz_t N) 
{
	initAll();
	fromAffine(x,y,N);
}

ExtendedPoint::ExtendedPoint(mpz_t N) 
{
	initAll();
	infinity(N);
}

void ExtendedPoint::infinity(mpz_t N)
{
	mpz_t x,y;
	mpz_init_set_ui(x,0);
	mpz_init_set_ui(y,1);

	fromAffine(x,y,N);
	mpz_clrs(x,y);
}

void ExtendedPoint::fromAffine(mpz_t x,mpz_t y,mpz_t N)
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

bool ExtendedPoint::toAffine(mpz_t x,mpz_t y,mpz_t N,mpz_t invB,mpz_t fact)
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

int readCurves(string file,mpz_t N,ExtendedPoint** pInit,bool& minus1)
{
	ifstream fp;
	if (file.length() < 2) 
	{
		file = "curves_twisted.txt";
		cout << "INFO: Defaulting to " << file  << endl;	
	}
	
	// Pokud se nepodari otevrit soubor, skonci
	fp.open(file);
	if (!fp.is_open())
	{
		cout << "ERROR: Cannot open file." << endl;
		return 0;
	} 
	
	// Inicializace racionalnich souradnic
	mpq_t qX,qY;
	mpq_intz(qX,qY);

	// Inicializace celociselnych souradnic
	mpz_t zX,zY;
	mpz_init_set_ui(zX,0);
	mpz_init_set_ui(zY,0);
	
	string ln;
	vector<ExtendedPoint> v;
	minus1 = true;
	cout << "Loading curves..." << endl;
	while (getline(fp,ln))
	{
		// Preskoc segment, ktery nezacina #
		if (ln.find("#") == string::npos) continue;
		
		// Je to prekroucena Edwardsova krivka s a = -1 ? 
		fp >> ln;
		if (v.size() > 0 && (ln == "-1") != minus1)
		{
			cout << "ERROR: Cannot mix curves with a = 1 and curves with a = -1." << endl; 
		
			fp.close();
			mpz_clrs(zX,zY);
			mpq_clrs(qX,qY);

			return 0;
		}
		minus1 = (ln == "-1");
		
		// Precti racionalni X-ovou souradnici a zkrat
		fp >> ln;
		mpq_set_str(qX,ln.c_str(),10);
		mpq_canonicalize(qX);
		
		// Precti racionalni Y-ovou souradnici a zkrat
		fp >> ln;
		mpq_set_str(qY,ln.c_str(),10);
		mpq_canonicalize(qY);

		// Pokus se X-ovou a Y-ovou souradnici rekudovat modulo N
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

		// Vytvor bod v Extended souradnicích z redukovanych afinnich bodu modulo N
		v.push_back(ExtendedPoint(zX,zY,N)); 
	}

	// Prekopiruj body z vektoru do pameti
	*pInit = new ExtendedPoint[v.size()];
	std::copy(v.begin(),v.end(),*pInit);

	// Vycisti pamet
	fp.close();
	mpz_clrs(zX,zY);
	mpq_clrs(qX,qY);

	return (int)v.size();
}
