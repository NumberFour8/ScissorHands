#include "generators.h"
#include "../helpers.h"

FileGenerator::FileGenerator(string filename)
 : Generator()
{
	fp.open(filename);
	S = 0;
}

FileGenerator::~FileGenerator()
{
	if (fp.is_open()) fp.close();
}

void FileGenerator::reset()
{
	if (fp.is_open()) 
	{
	  S = 0;
	  fp.clear();
	  fp.seekg(0,ios::beg);
	}
}

bool FileGenerator::next(ReducedPoint& P,const mpz_t zN)
{	
	string ln,t1;	

	if (!fp.is_open()) 
	{
		cout << "ERROR: Cannot open file." << endl;
		return false;
	}

	do { 
		if (!getline(fp,ln))
		  return false;  // Jsme na konci souboru
	} // Preskoc segment, ktery nezacina #
	while (ln.find("#") == string::npos);
	
	// Je to prekroucena Edwardsova krivka s a = -1 ? 
	fp >> ln;
	if (ln != "-1" && ln != "1")
	{
		cout << "ERROR: Unsupported curve type." << endl;
		return false; 
	}

	A = (ln == "-1" ? -1 : 1);

	// Precti racionalni souradnice 
	fp >> ln;
	fp >> t1;
	
	// Pokus se souradnice rekudovat modulo N
	P.reduce(RationalPoint(ln,t1),zN);

	return true; 
}
