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
	bool useDoubleAndAdd;
	bool verbose;
	bool noLCM;
	bool exitOnFinish;
	
	progArgs() : verbose(false), noLCM(false), useDoubleAndAdd(false) { }
};

// Precte parametry programu
void parseArguments(int argc,char** argv,progArgs& args)
{
	po::options_description desc("List of supported options");
	desc.add_options()
		("help", "Print usage information.")
		("verbose", po::value<bool>(&args.verbose)->default_value(false),
			"More verbose output.")
		("double-add", po::value<bool>(&args.useDoubleAndAdd)->default_value(false),
			"Use double-and-add instead sliding window.")
		("dont-compute-bound", po::value<bool>(&args.noLCM)->default_value(false),
			"True for s = B1, false for s = lcm(1,2...B1).")
		("no-restart", po::value<bool>(&args.exitOnFinish)->default_value(false),
			"When set, program terminates automatically after finishing.")
		("N-to-factor", po::value<string>(),
			"Number to factor.")
		("curve-file", po::value<string>(),
			"Path to file containing curves used for factoring.")
		("stage1-bound", po::value<unsigned int>(&args.B1)->default_value(4096),
			"Bound for ECM stage 1.")
		("window-size", po::value<unsigned short>(&args.windowSize)->default_value(4),
			"Size of sliding window or NAF width (in case of double-and-add).");
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);
	
	if (vm.count("N-to-factor"))
	  args.N = vm["N-to-factor"].as<string>();
	if (vm.count("curve-file"))
	  args.curveFile = vm["curve-file"].as<string>();
	if (vm.count("help"))
	  cout << endl << desc << endl << "-----------------------------------------" << endl;
}

// Zkontroluje úplnost parametrů, případně požádá o doplnění
void validateArguments(progArgs& args)
{
	const static boost::regex numVal("^[0-9]+$");
	if (args.N.empty() || !regex_match(args.N,numVal))
	{
		// Načíst N
		cout << "Enter N:" << endl;
		cin  >> args.N;
		cout << endl;
	}
	
	const fs::path p(curveFile);
	if (args.curveFile.empty() || !fs::exists(p) || !fs::is_regular_file(p))
	{
	  	// Načíst název souboru s křivkami
		cout << "Enter path to curve file:" << endl;
		cin  >> args.curveFile;
		cout << endl;
	}
	if (args.B1 <= 2)
	{
		// Načíst hranici první fáze
		cout << "Enter stage 1 bound B1:" << endl;
		cin  >> args.B1;
		cout << endl;
	}
	if (args.windowSize < 2)
	{
		// Velikost okna
		cout << "Enter window size:" << endl;
		cin  >> args.windowSize;
		cout << endl;
	}
	validateArguments(args);
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
	set<factor,bool(*)(factor x, factor y)> foundFactors([](factor x, factor y){ return x.fac < y.fac; });

	// Jsou-li předány parametry, použij je. Jinak se na ně zeptej.
	parseArguments(argc,argv,args);
	validateArguments(args);

	// Inicializace N
	mpz_t zN;
	mpz_init_set_str(zN,args.N.c_str(),10);
	
	// Pokud je N prvočíslo, není co faktorizovat
	if (mpz_probab_prime_p(zF,25) != 0)
	{
		cout << "ERROR: N is almost surely prime." << endl;
		exitCode = 1;
		goto end;
	}

	// Inicializace proměnných
	ExtendedPoint infty(zN); 	// Neutrální prvek
	ComputeConfig ax(zN);    	// Pomocná struktura
	NAF S;					 	// NAF rozvoj
	mpz_t zS,zInvW,zX,zY,zF; 	// Pomocné proměnné
	cudaError_t cudaStatus;	 	// Proměnná pro chybové kódy GPU
	ExtendedPoint *PP;		 	// Adresa všech bodů
	bool minusOne;			 	// Pracujeme s křivkami s a =-1 ?

	restart_bound:

	PP = NULL;
	read_curves = readCurves(args.curveFile,zN,&PP,minusOne);

	// Zkontroluj počet načtených křivek
	if (read_curves <= 0)
	{
		cout << "ERROR: No curves read." << endl;
		exitCode = 1;
		goto end;
	}
	cout << "Loaded " << read_curves << " curves with a = " << (minusOne ? "-1" : "1")  << endl << endl;


	// Spočti S = lcm(1,2,3...,B1) a jeho NAF rozvoj
	mpz_init(zS);
	cout << "Computing coefficient..." << endl;
	
	if (args.noLCM)
	  mpz_set_ui(zS,args.B1);
	else lcmToN(zS,args.B1);
	
	S.initialize(zS,(unsigned char)args.windowSize);
	mpz_clear(zS);	
	
	cout << endl << "Trying to factor " << args.N << " with B1 = "<< args.B1 << " using " << read_curves << " curves..." << endl << endl;

	// Nastavit hodnoty do konfigurace
	ax.windowSz  = args.windowSize;
	ax.nafLen    = S.l;
	ax.numCurves = read_curves;
	ax.minus1	 = minusOne; 
	ax.useDblAdd = args.useDoubleAndAdd;

	// Proveď výpočet
	cudaStatus = compute(ax,&infty,PP,S);
    if (cudaStatus != cudaSuccess) 
    {
        cout << "ERROR: CUDA compute failed!" << endl;
		exitCode = 1;
        goto end;
    }
	cout << "Computation finished." << endl;

	// Inicializace pomocných proměnných
	mpz_intz(zInvW,zX,zY,zF);
	
	// Spočti 2^(-W) mod N 
	mpz_ui_pow_ui(zInvW, 2, SIZE_DIGIT); 
	mpz_invert(zInvW, zInvW, zN);

	// Analyzuj výsledky
	foundFactors.clear();
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
		   string fact = mpz_to_string(zF);
		   cout << "Factor found: " << fact << endl;
		   foundFactors.insert(factor(fact,mpz_probab_prime_p(zF,25) != 0,i));
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
	else cout << endl << "No factors found." << endl << endl;

	// Vyčisti proměnné
	mpz_clrs(zInvW,zX,zY);
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
