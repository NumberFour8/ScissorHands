#include <stdio.h>
#include <iomanip>
#include <sstream> 
#include <set>

#include <boost/regex.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

#include "kernel.h"

const unsigned int MAX_STAGE1_BOUND = 500000;

// Struktura uchovavajici konfiguraci prectenou z parametru nebo ze vstupu
struct progArgs {
	string N;
	vector<string> curveFiles;
	vector<unsigned int> B1;
	unsigned int curB1,curCur;
	unsigned short windowSize;
	bool verbose;
	bool noLCM;
	bool exitOnFinish;
	bool greedy;
	unsigned short whichDevice;
	string outputFile;
	
	progArgs() 
	: verbose(false), noLCM(false), greedy(false), exitOnFinish(false), curB1(0), curCur(0), windowSize(0), whichDevice(0)
	{ }
};

// Přečte parametry programu
void parseArguments(int argc,char** argv,progArgs& args)
{
	po::options_description desc("List of supported options");
	desc.add_options()
		("help,h", "Print usage information.")
		("verbose,v", "More verbose output.")
		("dont-compute-bound,x", "Coefficient s = B1 when set, otherwise s = lcm(1,2...B1).")
		("no-restart,e", "When set, program terminates automatically after finishing.")
		("device-id,D", po::value<unsigned short>(&args.whichDevice)->default_value(0),
			"ID of a CUDA device used for computation.")
		("N-to-factor,N", po::value<string>(),
			"Number to factor.")
		("curve-files,f", po::value<vector<string>>()->multitoken(),
			"Path to file(s) containing curves used for factoring.")
		("stage1-bound,B", po::value<vector<unsigned int>>()->multitoken(),
			"Bound for ECM stage 1 or range set by 'start stride end'")
		("greedy,g", "If set, wait for input of another N to factor after finishing.") 
		("window-size,W", po::value<unsigned short>(&args.windowSize)->default_value(4),
			"Size of sliding window.")
		("output-file,o",po::value<string>()->default_value("factors-of"),
			"File name where to output all found prime factors.");
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);
	
	args.verbose		 = vm.count("verbose") != 0;
	args.noLCM			 = vm.count("dont-compute-bound") != 0;
	args.exitOnFinish	 = vm.count("no-restart") != 0 || vm.count("greedy") != 0;;
	args.greedy			 = vm.count("greedy") != 0;

	if (vm.count("N-to-factor"))
	  args.N = vm["N-to-factor"].as<string>();
	if (vm.count("curve-files"))
	  args.curveFiles = vm["curve-files"].as<vector<string>>();
	if (vm.count("stage1-bound"))
	  args.B1 = vm["stage1-bound"].as<vector<unsigned int>>();
	if (vm.count("output-file"))
	  args.outputFile = vm["output-file"].as<string>();
	if (vm.count("help"))
	  cout << endl << desc << endl << "-----------------------------------------" << endl << endl;
}

// Zkontroluje úplnost parametrů, případně požádá o doplnění
int validateArguments(progArgs& args)
{
	bool recheck = false; 
	const static boost::regex numVal("^[0-9]+$");
	if (args.N.empty() || !regex_match(args.N,numVal))
	{
		// Načíst N
		cout << "Enter N:" << endl;
		if (!getline(cin,args.N))
			return 1;
		cout << endl;
		recheck = true; 
	}
	
	if (args.curveFiles.size() > 0) 
	{
		for (unsigned int i = 0;i < args.curveFiles.size();++i)
		{
		   // Zkontrolovat a nahradit krátké zástupce souborů
		   if (args.curveFiles[i].length() <= 3)
		   {
			  args.curveFiles[i] =  args.curveFiles[i] == "E" ? "curves_edwards.txt" :
									args.curveFiles[i] == "M" ? "curves_mixed.txt"   : "curves_twisted.txt";
			  cout << "INFO: Defaulting to " << args.curveFiles[i]  << endl;	
		   }
	
		   // Zkontrolovat existenci souborů
		   fs::path p(args.curveFiles[i]);
		   if (!fs::exists(p) || !fs::is_regular_file(p))
		   {
			  // Načíst znovu název souboru s křivkami, pokud je neplatný
			  cout << "Path to curve file " << args.curveFiles[i] << " is invalid, please re-enter path:" << endl;
			  cin  >> args.curveFiles[i];
			  cout << endl;
			  recheck = true;
		   }
		}
	}
	else 
	{ 
		// Načíst soubory s křivkami, není-li žádný uveden
		string files;
		cout << "Enter path to at least one curve file (comma separated):" << endl;
		cin  >>  files;
		cout << endl;
		
		boost::split(args.curveFiles,files,boost::is_any_of(";,"));
		recheck = true;
	}

	if (args.B1.size() == 1 || args.B1.size() == 3)
	{
		// Zkontrolovat velikost první hodnoty
		if (args.B1[0] < 2 || args.B1[0] > MAX_STAGE1_BOUND)
		{
			cout << "Enter stage 1 bound B1:" << endl;
			cin  >> args.B1[0];
			cout << endl;
			recheck = true;
		}
		if (args.B1.size() == 3)
		{
			// Zkontrolovat velikost poslední hodnoty
			if (args.B1[0] >= args.B1[2] || args.B1[2] > MAX_STAGE1_BOUND) 
			{
				cout << "Enter stage 1 bound end:" << endl;
				cin  >> args.B1[2];
				cout << endl;
				recheck = true;
			} // Zkontrolovat velikost kroku
			else if (args.B1[1] > args.B1[2] || args.B1[1] < 128)
			{
				cout << "Enter stage 1 bound stride:" << endl;
				cin  >> args.B1[1];
				cout << endl;
				recheck = true;
			}
		}
	}
	else 
	{
		// Načíst hranice první fáze
		string bounds;
		cout << "Enter stage 1 bound B1 or Start,Stride,End for range:" << endl;
		cin  >> bounds;
		cout << endl;
		
		vector<string> bs;
		boost::split(bs,bounds,boost::is_any_of(";,:"));
		
		std::transform(bs.begin(),bs.end(),std::back_inserter(args.B1),
					   [](const string& s) { return std::stoi(s); });
		
		recheck = true;
	}

	if (args.windowSize < 2 || args.windowSize > 15)
	{
		// Velikost okna
		cout << "Enter window size:" << endl;
		cin  >> args.windowSize;
		cout << endl;
		recheck = true;
	}
	
	args.curB1  = args.B1[0];
	args.curCur = 0;
	if (recheck) return validateArguments(args);
	return 0;
}

// Uloží a přepíše výstupní soubor s nalezenými faktory
void savePrimeFactors(string fileName,int id,stringstream& primeStream)
{
	string fname = (boost::format("%s-%d.txt") % fileName % id).str();
	ofstream pr(fname,ofstream::out | ofstream::trunc);
	pr << primeStream.str();
	pr.close();
	cout << "All found prime factors have been written to file: " << fname << endl;
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
	cout << "ScissorHands ECM Factorization" << endl;
	cout << "by Lukas Pohanka, FNSPE CTU, 2014" << endl;
	cout << endl;
	
	progArgs args;
	stringstream primeStream;
	int exitCode = 0,lastB1 = 0,runNum = 0,factorCount = 0,Ncount = 1;
	float cudaTimeCounter = 0;
	bool useMixedStrategy = false,fullFactorizationFound = false;
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
	
	primeStream << "# Found prime factors of " << args.N << "\n";
	primeStream << "# FOUND FACTOR, CURVE ID\n";
	
	// Inicializace proměnných
	ExtendedPoint infty;			 // Neutrální prvek
	ComputeConfig ax;	    		 // Pomocná struktura
	NAF S;					 		 // NAF rozvoj
	mpz_t zS,zInvW,zX,zY,zF,zChk; 	 // Pomocné proměnné
	cudaError_t cudaStatus;	 		 // Proměnná pro chybové kódy GPU
	ExtendedPoint *PP;		 		 // Adresa všech bodů
	int read_curves,edwards,twisted; // Počty načtených typů křivek
	int bitCountN;					 // Počet bitů N
	computeStrategy strategy;		 // Strategie výpočtu

	restart_bound:
	
	// Pokud je N prvočíslo, není co faktorizovat
	if (is_almost_surely_prime(zN) || mpz_cmp_ui(zN,1) == 0 || mpz_cmp_ui(zN,0) == 0)
	{
		cout << "ERROR: N equals 0,1 or is almost surely a prime." << endl;
		exitCode = 1;
		goto end;
	}
	
	// Kontrola velikosti N
	bitCountN = (int)mpz_sizeinbase(zN,2);
	if (bitCountN > NB_DIGITS*SIZE_DIGIT)
	{
		cout << "ERROR: Cannot factor numbers longer than " << NB_DIGITS*SIZE_DIGIT << " bits." << endl;
		exitCode = 1;
		goto end;
	}
	cout << "NOTE: N is " << bitCountN << " bits long." << endl;

	infty.infinity(zN);
	ax.initialize(zN);

	PP = NULL;
	read_curves = edwards = twisted = 0;
	strategy	= computeStrategy::csNone;
	
	// Načti křivky a zvol vhodnou strategii
	strategy	= readCurves(args.curveFiles[args.curCur],zN,&PP,edwards,twisted,read_curves);
	if (strategy == computeStrategy::csNone)
	{
		cout << "ERROR: No suitable compute strategy found." << endl;
		exitCode = 1;
		goto end;
	} 
	
	// Spočti S = lcm(1,2,3...,B1) a jeho NAF rozvoj, pokud se B1 změnilo
	if (lastB1 != args.curB1)
	{
		mpz_init(zS);
		cout << "Recomputing coefficient..." << endl;
		
		if (args.noLCM)
		{
		  cout << "NOTE: Using bound B1 as a coefficient directly." << endl; 
		  mpz_set_ui(zS,args.curB1);
		}
		else lcmToN(zS,args.curB1);
		lastB1 = args.curB1;

		S.initialize(zS,2);
		mpz_clear(zS);	
	}
	else cout << "NOTE: B1 hasn't changed." << endl;  
	
	cout << endl << "Trying to factor " << mpz_to_string(zN) << " with B1 = "<< args.curB1 << " using " << read_curves << " curves..." << endl << endl;

	// Nastavit hodnoty do konfigurace
	ax.windowSz  = args.windowSize;
	ax.nafLen    = S.l;
	ax.numCurves = read_curves;
	ax.minus1	 = (strategy == computeStrategy::csTwisted);
	ax.deviceId	 = args.whichDevice;
	ax.cudaRunTime = 0;
	runNum++;

	// Proveď výpočet
	if (strategy == computeStrategy::csMixed)
	{
		cout << "NOTE: Using mixed compute strategy." << endl;
		cudaStatus = computeMixed(ax,&infty,PP,S);
	}
	else 
	{
		cout << "NOTE: Using single compute strategy." << endl;
		cudaStatus = computeSingle(ax,&infty,PP,S);
	}
	
	if (cudaStatus != cudaSuccess) 
    {
        cout << "ERROR: CUDA compute failed!" << endl;
		exitCode = 1;
		delete[] PP;
        goto end;
    }
	cout << "Computation finished." << endl << endl;
	cudaTimeCounter += ax.cudaRunTime;

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
		if (args.verbose) cout << "Curve #" << i+1 << ":\t"; 
		if (PP[i].toAffine(zX,zY,zN,zInvW,zF)) 
		{
			if (args.verbose)
			{
			  cout << "No factor found." << endl;
			  cout << endl << "sP = (" << mpz_to_string(zX) << "," << mpz_to_string(zY) << ")" << endl;
			}
		}
		else if (mpz_cmp_ui(zF,0) != 0) // Máme faktor!
		{
		   bool isPrime = is_almost_surely_prime(zF);
		   string fact  = mpz_to_string(zF);

		   if (args.verbose) cout << "Factor found: " << fact << endl;
		   if (foundFactors.insert(factor(fact,isPrime,i)).second && isPrime) 
		   {
			 mpz_mul(zChk,zChk,zF);
		     primeStream << fact << ", " << i;
			 if (strategy == computeStrategy::csEdwards)
				primeStream << "E\n";
			 else if (strategy == computeStrategy::csTwisted)
				primeStream << "T\n";
			 else primeStream << "M\n";
			 factorCount++;
		   }
		}
		else if (args.verbose) cout << "Error during conversion." << endl;
		if (args.verbose) cout << endl << "------------" << endl;
    }
	
	// Vypiš všechny nalezené faktory
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
			fullFactorizationFound = true;
		}
		else {
			cout << "REMAINING UNFACTORED PART: " << mpz_to_string(zChk); 
			cout << " (" << mpz_sizeinbase(zChk,2) << " bits)" << endl;
			mpz_set(zN,zChk);
		}
	}
	else 
	{
		cout << mpz_to_string(zN) << " HAS BEEN FACTORED TOTALLY!" << endl;
		args.exitOnFinish = true; 
		fullFactorizationFound = true;
	}

	// Vyčisti proměnné
	mpz_clrs(zInvW,zX,zY,zChk);
	delete[] PP;

	// Kontrola,zda ještě můžeme použít nějaký soubor s křivkami
	if (args.curCur < args.curveFiles.size()-1 && !fullFactorizationFound)
	{
		args.curCur += 1;
		cout << endl << "NOTE: Trying a different curve file." << endl;
		goto restart_bound;
	}

	// Kontrola, zda máme restartovat bez zeptání a zvýšit B1
	if (args.B1.size() > 1 && args.curB1 <= args.B1[2]-args.B1[1] && !fullFactorizationFound)
	{
		args.curB1 += args.B1[1];
		args.curCur = 0;
		cout << "NOTE: B1 has been automatically incremented to: " << args.curB1 << endl;
		goto restart_bound;
	}

	// Finální menu
	while (!args.exitOnFinish)
	{
	   cout << endl << "Type : " << endl;
	   cout << "'r' to restart with new configuration." << endl;
	   cout << "'p' to print all found prime factors so far." << endl;
	   cout << "'a [number]' to increase B1 by [number] and restart with the same curve file." << endl; 
	   cout << "or 'q' to quit." << endl << endl;
	   cin  >> c;
	   cout << endl;
	   if (c == 'q') break;
	   if (c == 'a')
	   {
		  unsigned int B1_inc;
		  cin >> B1_inc;
		  
		  args.B1[0] += B1_inc;
		  validateArguments(args);

		  cout << "B1 has been incremented to: " << args.curB1 << endl;
		  goto restart_bound;
	   }
	   else if (c == 'p')
	   {
		  cout << endl;
		  cout << primeStream.str() << endl << endl;
	   }
	   else if (c == 'r')
	   {
		  args.B1.clear();
		  args.curveFiles.clear();

		  validateArguments(args);
		  goto restart_bound;
	   }
	}
	
	// Ulož výstup a vypiš celkový čas běhu
	primeStream << "Found prime factors: " << factorCount;
	savePrimeFactors(args.outputFile,Ncount,primeStream);
	cout << "Total GPU running time is : " << setprecision(3) << (cudaTimeCounter/60000) << " minutes." << endl;

	// Jsme-li v hladovém módu, chtěj další číslo k faktorizaci
	if (args.greedy)
	{
		args.N = "";
		lastB1 = runNum = factorCount = 0;
		cudaTimeCounter = 0;
		primeStream.str(string(""));
		fullFactorizationFound = false; 
		
		cout << endl; 
		if (validateArguments(args) != 0) goto end;
		Ncount++;
		goto restart_bound;
	}

	end:
	mpz_clear(zN);

	return exitCode;
}
