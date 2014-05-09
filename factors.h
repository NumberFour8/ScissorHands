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

typedef bool(*factor_comparator)(factor, factor);

class FoundFactors {
private:
	
	// Množina nalezených faktorů s vlastním uspořádnáním
	set<factor,factor_comparator> newFoundFactors;
	int totalPrimesFound,currentSessionId,primesFoundInSession;
	
	stringstream primeStream;
	string fact;
	mpz_t zChk;
	
public:

	FoundFactors() : totalPrimesFound(0), currentSessionId(0), primesFoundInSession(0),
		newFoundFactors([](factor x, factor y){ return x.fac < y.fac; })
	{
		mpz_init(zChk);
	}
	
	~FoundFactors()
	{
		mpz_clear(zChk);
	}
	
	void startNewSession(string zN)
	{
		currentSessionId++;
		
		primeStream.str(string(""));
		primeStream << "# Session ID: " << currentSessionId << "\n";
		primeStream << "# Found prime factors of " << zN << "\n";
		primeStream << "# FOUND FACTOR, CURVE ID\n";
		
		mpz_set_ui(zChk,1);
		newFoundFactors.clear();
		primesFoundInSession = 0;
	}

	string getLastFactor()
	{
		return fact;
	}
	
	int primesFoundInAllSessions()
	{
		return totalPrimesFound;
	}
	
	int primesFoundInCurrentSession()
	{
		return primesFoundInSession;
	}

	// Zkontroluj a přidej nový faktor
	bool handlePotentialFactor(mpz_t zF,int curveId)
	{
		fact.clear();
		if (mpz_cmp_ui(zF,0) != 0) // Máme faktor!
		{
		   bool isPrime = is_almost_surely_prime(zF);
		   fact  = mpz_to_string(zF);

		   if (newFoundFactors.insert(factor(fact,isPrime,curveId)).second && isPrime) 
		   {
			 mpz_mul(zChk,zChk,zF);
			 primeStream << fact << ", " << curveId << "\n";
			 
			 totalPrimesFound++;
			 primesFoundInSession++;
		   }
		   return true;
		}
		else return false;
	}
	
	// Vypiš všechny nově nalezené faktory
	void printNewFoundFactors(int runNum)
	{
		if (newFoundFactors.size() > 0)
		{
		  cout << endl << newFoundFactors.size() << " FACTORS FOUND IN RUN #" << runNum << ":" << endl << endl;
		  std::for_each(newFoundFactors.begin(),newFoundFactors.end(),
			[](const factor f)
			{ 
				cout << (f.prime ? "Prime:\t\t" : "Composite:\t") << f.fac << " (#" << f.curveId << ")" <<  endl; 
			}
		   );
		}
		else cout << endl << "NO FACTORS FOUND." << endl << endl;
	}
	
	void printPrimesFromSession()
	{
		cout << primeStream.str() << endl;
	}
	
	// Uloží a přepíše výstupní soubor s nalezenými faktory
	void savePrimesFromSession(string fileName,bool append)
	{
		ios_base::openmode md = ofstream::out | ofstream::trunc;
		string   fn	= (boost::format("%s-%d.txt") % fileName % currentSessionId).str();
		if (append)
		{
		  md = ofstream::out | ofstream::app;
		  fn = (boost::format("%s.txt") % fileName).str();
		}
			
		ofstream pr(fn,md);
		primeStream << boost::format("\n# Found prime factors: %d\n# -------------------------------------\n\n") % primesFoundInSession;
		pr << primeStream.str();
		primeStream.str(string(""));
		
		pr.close();
		cout << "All found prime factors have been written to file: " << fn << endl;
	}
	
	bool tryFinishSession(mpz_t zN)
	{
		newFoundFactors.clear();
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
				mpz_set_ui(zChk,1);
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
