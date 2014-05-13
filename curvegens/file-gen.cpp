#include "generators.h"
#include "../helpers.h"

FileGenerator::FileGenerator(string filename,unsigned int cacheSize)
 : Generator(cacheSize)
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

void FileGenerator::revert()
{
	Generator::revert();	
}

bool FileGenerator::next(RationalPoint& P)
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

	P.minus1 = (ln == "-1");

	// Precti racionalni souradnice 
	fp >> ln;
	fp >> t1;

	P.set(ln,t1);
	
	return true; 
}
