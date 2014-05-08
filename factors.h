#ifndef FOUNDFACTORS_H
#define FOUNDFACTORS_H

#include <sstream> 
#include <set>

#include <boost/format.hpp>

#include "helpers.h"

// Struktura obsahující informace o získaném faktoru
struct factor {
		string fac;
		bool prime;
		unsigned int curveId;
		factor(string f,bool p,unsigned int c) : fac(f),prime(p),curveId(c) {}
};


class FoundFactors {
private:
	
	// Množina nalezených faktorů s vlastním uspořádnáním
	typedef bool(*factor_comparator)(factor, factor);
	set<factor,factor_comparator> foundFactors([](factor x, factor y){ return x.fac < y.fac; });
	int totalFound;
	
	stringstream primeStream;
	string fact;
	mpz_t zChk;
	
public:


	FoundFactors(string initialN) : totalFound(0)
	{
		mpz_init(zChk);
		primeStream << "# Found prime factors of " << initialN << "\n";
		primeStream << "# FOUND FACTOR, CURVE ID\n";
	}
	
	~FoundFactors()
	{
		mpz_clear(zChk);
	}

	string getLastFactor()
	{
		return fact;
	}

	void clearSession()
	{
		foundFactors.clear();
		mpz_set_ui(zChk,1);
	}

	// Zkontroluj a přidej nový faktor
	bool handlePotentialFactor(mpz_t zF,int curveId)
	{
		fact.clear();
		if (mpz_cmp_ui(zF,0) != 0) // Máme faktor!
		{
		   bool isPrime = is_almost_surely_prime(zF);
		   fact  = mpz_to_string(zF);

		   if (foundFactors.insert(factor(fact,isPrime,curveId)).second && isPrime) 
		   {
			 mpz_mul(zChk,zChk,zF);
			 primeStream << fact << ", " << curveId << "\n";
			 totalFound++;
		   }
		   return true;
		}
		else return false;
	}
	
	// Vypiš všechny nalezené faktory
	void printAllFactors()
	{
		if (foundFactors.size() > 0)
		{
		  cout << endl << foundFactors.size() << " FACTORS FOUND IN RUN #" << runNum << ":" << endl << endl;
		  std::for_each(foundFactors.begin(),foundFactors.end(),
			[](const factor f)
			{ 
				cout << (f.prime ? "Prime:\t\t" : "Composite:\t") << f.fac << " (#" << f.curveId << ")" <<  endl; 
			}
		   );
		}
		else cout << endl << "NO FACTORS FOUND." << endl << endl;
	}
	
	// Uloží a přepíše výstupní soubor s nalezenými faktory
	void savePrimeFactors(string fileName,int fileId,bool append)
	{
		ios_base::openmode md = ofstream::out;
		string fname = "";
		if (append)
		{
			fname = (boost::format("%s-all.txt") % fileName).str();
			md |= ofstream::app;
		}
		else 
		{
			fname = (boost::format("%s-#%d.txt") % fileName % fileId).str();
			md |= ofstream::trunc;
		}
		
		ofstream pr(fname,md);
		
		primeStream << boost::format("\n# Found prime factors: %d\n# -------------------------------------\n\n") % totalFound;
		pr << primeStream.str();
		pr.close();
		cout << "All found prime factors have been written to file: " << fname << endl;
	}
	
	bool reduceComposite(mpz_t zN)
	{
		if (mpz_cmp(zChk,zN) != 0)
		{
			mpz_divexact(zChk,zN,zChk);
			if (is_almost_surely_prime(zChk))
			{
				cout << "REMAINING UNFACTORED PART " << mpz_to_string(zChk) << " IS A PRIME." << endl;
				cout << mpz_to_string(zN) << " HAS BEEN FACTORED TOTALLY!" << endl;
				return true;
			}
			else {
				cout << "REMAINING UNFACTORED PART: " << mpz_to_string(zChk); 
				cout << " (" << mpz_sizeinbase(zChk,2) << " bits)" << endl;
				mpz_set(zN,zChk);
				return false;
			}
		}
		else 
		{
			cout << mpz_to_string(zN) << " HAS BEEN FACTORED TOTALLY!" << endl;
			return true;
		}
	}
	
};


#endif
