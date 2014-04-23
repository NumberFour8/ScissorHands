#include <stdio.h>
#include <set>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "kernel.h"

// Struktura uchovavajici konfiguraci prectenou z parametru nebo ze vstupu
struct progArgs {
	string N;
	string curveFile;
	unsigned int B1;
	unsigned short windowSize;
	bool verbose;
	bool noLCM;
	bool exitOnFinish;
	
	progArgs() 
	 : verbose(false), noLCM(false), B1(0), windowSize(0) 
	{ }
};

// Precte parametry programu
void parseArguments(int argc,char** argv,progArgs& args)
{
	po::options_description desc("List of supported options");
	desc.add_options()
		("help", "Print usage information.")
		("verbose", "More verbose output.")
		("dont-compute-bound", "Coefficient s = B1 when set, otherwise s = lcm(1,2...B1).")
		("no-restart", "When set, program terminates automatically after finishing.")
		("N-to-factor", po::value<string>(),
			"Number to factor.")
		("curve-file", po::value<string>(),
			"Path to file containing curves used for factoring.")
		("stage1-bound", po::value<unsigned int>(&args.B1),
			"Bound for ECM stage 1.")
		("window-size", po::value<unsigned short>(&args.windowSize),
			"Size of sliding window (or NAF width in case of double-and-add).");
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);
	
	args.verbose		 = vm.count("verbose") != 0;
	args.noLCM			 = vm.count("dont-compute-bound") != 0;
	args.exitOnFinish	 = vm.count("no-restart") != 0;

	if (vm.count("N-to-factor"))
	  args.N = vm["N-to-factor"].as<string>();
	if (vm.count("curve-file"))
	  args.curveFile = vm["curve-file"].as<string>();
	if (vm.count("help"))
	  cout << endl << desc << endl << "-----------------------------------------" << endl << endl;
}

// Zkontroluje úplnost parametrů, případně požádá o doplnění
void validateArguments(progArgs& args)
{
	bool recheck = false; 
	const static boost::regex numVal("^[0-9]+$");
	if (args.N.empty() || !regex_match(args.N,numVal))
	{
		// Načíst N
		cout << "Enter N:" << endl;
		cin  >> args.N;
		cout << endl;
		recheck = true; 
	}
	
	if (!args.curveFile.empty() && args.curveFile.length() < 2) 
	{
		args.curveFile = args.curveFile == "E" ? "curves_edwards.txt" :
					     args.curveFile == "M" ? "curves_mixed.txt"   : "curves_twisted.txt";
		cout << "INFO: Defaulting to " << args.curveFile  << endl << endl;	
	}
	else if (args.curveFile.empty())
	{ 
		cout << "Enter path to curve file:" << endl;
		cin  >> args.curveFile;
		cout << endl;
		recheck = true;
	}
	else 
	{ 
		const fs::path p(args.curveFile);
		if (!fs::exists(p) || !fs::is_regular_file(p))
		{
	  		// Načíst název souboru s křivkami
			cout << "Enter path to curve file:" << endl;
			cin  >> args.curveFile;
			cout << endl;
			recheck = true;
		}
	}

	if (args.B1 <= 2)
	{
		// Načíst hranici první fáze
		cout << "Enter stage 1 bound B1:" << endl;
		cin  >> args.B1;
		cout << endl;
		recheck = true;
	}

	if (args.windowSize < 2)
	{
		// Velikost okna
		cout << "Enter window size:" << endl;
		cin  >> args.windowSize;
		cout << endl;
		recheck = true;
	}
	if (recheck) validateArguments(args);
}

// Struktura obsahující informace o získaném faktoru
struct factor {
		string fac;
		bool prime;
		unsigned int curveId;
		factor(string f,bool p,unsigned int c) : fac(f),prime(p),curveId(c) {}
};

int main(int argc,char** argv)
{
	cout << "ECM using Twisted Edwards curves" << endl;
	
	progArgs args;
	int read_curves = 0,exitCode = 0;
	char c = 0;

	// Množina nalezených faktorů s vlastním uspořádnáním
	typedef bool(*factor_comparator)(factor, factor);
	set<factor,factor_comparator> foundFactors([](factor x, factor y){ return x.fac < y.fac; });

	// Jsou-li předány parametry, použij je. Jinak se na ně zeptej.
	parseArguments(argc,argv,args);
	validateArguments(args);

	// Inicializace N
	mpz_t zN;
	mpz_init_set_str(zN,args.N.c_str(),10);
	
	// Inicializace proměnných
	ExtendedPoint infty;			// Neutrální prvek
	ComputeConfig ax;	    		// Pomocná struktura
	NAF S;					 		// NAF rozvoj
	mpz_t zS,zInvW,zX,zY,zF,zChk; 	// Pomocné proměnné
	cudaError_t cudaStatus;	 		// Proměnná pro chybové kódy GPU
	ExtendedPoint *PP;		 		// Adresa všech bodů

	restart_bound:
	
	// Pokud je N prvočíslo, není co faktorizovat
	if (is_almost_surely_prime(zN) || mpz_cmp_ui(zN,1) == 0 || mpz_cmp_ui(zN,0) == 0)
	{
		cout << "ERROR: N equals 0,1 or is almost surely a prime." << endl;
		exitCode = 1;
		goto end;
	}

	infty.infinity(zN);
	ax.initialize(zN);

	PP = NULL;
	read_curves = readCurves(args.curveFile,zN,&PP);

	// Zkontroluj počet načtených křivek
	if (read_curves <= 0)
	{
		cout << "ERROR: No curves read." << endl;
		exitCode = 1;
		goto end;
	}
	else if (read_curves < CURVES_PER_BLOCK*2)
	{
		cout << "ERROR: Minimum number of curves is " << CURVES_PER_BLOCK*2 << endl;
		exitCode = 1;
		goto end;
	}
	cout << "Loaded " << read_curves << " curves." << endl << endl;


	// Spočti S = lcm(1,2,3...,B1) a jeho NAF rozvoj
	mpz_init(zS);
	cout << "Computing coefficient..." << endl;
	
	if (args.noLCM)
	{
	  cout << "NOTE: Using bound B1 as a coefficient directly." << endl; 
	  mpz_set_ui(zS,args.B1);
	}
	else lcmToN(zS,args.B1);
	
	S.initialize(zS,2);
	mpz_clear(zS);	
	
	cout << endl << "Trying to factor " << args.N << " with B1 = "<< args.B1 << " using " << read_curves << " curves..." << endl << endl;

	// Nastavit hodnoty do konfigurace
	ax.windowSz  = args.windowSize;
	ax.nafLen    = S.l;
	ax.numCurves = read_curves;
	
	// Proveď výpočet
	cudaStatus = compute(ax,&infty,PP,S);
    if (cudaStatus != cudaSuccess) 
    {
        cout << "ERROR: CUDA compute failed!" << endl;
		exitCode = 1;
        goto end;
    }
	cout << "Computation finished." << endl << endl;

	// Inicializace pomocných proměnných
	mpz_intz(zInvW,zX,zY,zF);
	
	// Spočti 2^(-W) mod N 
	mpz_ui_pow_ui(zInvW, 2, SIZE_DIGIT); 
	mpz_invert(zInvW, zInvW, zN);

	// Analyzuj výsledky
	foundFactors.clear();
	mpz_init_set_ui(zChk,1);
	for (int i = 0; i < read_curves;++i)
	{
		cout << "Curve #" << i+1 << ":\t"; 
		if (PP[i].toAffine(zX,zY,zN,zInvW,zF)) 
		{
			cout << "No factor found." << endl;
			if (args.verbose)
			  cout << endl << "sP = (" << mpz_to_string(zX) << "," << mpz_to_string(zY) << ")" << endl;
		}
		else if (mpz_cmp_ui(zF,0) != 0) // Máme faktor!
		{
		   bool isPrime = is_almost_surely_prime(zF);
		   string fact  = mpz_to_string(zF);

		   cout << "Factor found: " << fact << endl;
		   if (foundFactors.insert(factor(fact,isPrime,i)).second && isPrime) 
			 mpz_mul(zChk,zChk,zF);
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
		{ 
			cout << (f.prime ? "Prime:\t" : "Composite:\t") << f.fac << " (#" << f.curveId << ")" <<  endl; 
		}
	   );
	}
	else cout << endl << "NO FACTORS FOUND." << endl << endl;

	// Poděl číslo N všemi nalezenými prvočíselnými faktory
	cout << endl; 
	if (mpz_cmp(zChk,zN) != 0)
	{
		mpz_divexact(zChk,zN,zChk);
		if (is_almost_surely_prime(zChk))
		{
			cout << "REMAINING UNFACTORED PART " << mpz_to_string(zChk) << " IS A PRIME." << endl;
			cout << mpz_to_string(zN) << " HAS BEEN FACTORED TOTALLY!" << endl;
			args.exitOnFinish = true; 
		}
		else {
			cout << "REMAINING UNFACTORED PART: " << mpz_to_string(zChk) << endl; 
			mpz_set(zN,zChk);
		}
	}
	else 
	{
		cout << mpz_to_string(zN) << " HAS BEEN FACTORED TOTALLY!" << endl;
		args.exitOnFinish = true; 
	}

	// Vyčisti proměnné
	mpz_clrs(zInvW,zX,zY,zChk);
	delete[] PP;

	while (!args.exitOnFinish)
	{
	   cout << endl << "Type 'r' to restart with new configuration or 'q' to quit..." << endl;
	   cin  >> c;
	   if (c == 'q') break;
	   else if (c == 'r')
	   {
		  args.B1 = args.windowSize = 0;
		  args.curveFile.clear();
		  
		  validateArguments(args);
		  goto restart_bound;
		}
	}
	
	end:
	mpz_clear(zN);

	return exitCode;
}
