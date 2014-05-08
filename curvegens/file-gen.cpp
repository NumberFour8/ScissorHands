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
	string ln;	
	bool ret = false;

	mpq_t qX,qY;
	mpq_intz(qX,qY);

	if (!fp.is_open()) 
	{
		cout << "ERROR: Cannot read from file." << endl;
		goto read_finish;
	}

	do { 
		if (!getline(fp,ln))
		  goto read_finish; // Jsme na konci souboru
	} // Preskoc segment, ktery nezacina #
	while (ln.find("#") == string::npos);
	
	// Je to prekroucena Edwardsova krivka s a = -1 ? 
	fp >> ln;
	if (ln != "-1" && ln != "1")
	{
		cout << "ERROR: Unsupported curve type." << endl;
		goto read_finish;
	}

	A = (ln == "-1" ? -1 : 1);

	// Precti racionalni X-ovou souradnici 
	fp >> ln;
	mpq_set_str(qX,ln.c_str(),10);
	
	// Precti racionalni Y-ovou souradnici 
	fp >> ln;
	mpq_set_str(qY,ln.c_str(),10);

	// Pokus se X-ovou a Y-ovou souradnici rekudovat modulo N
	reduce_rational_point(P.X,P.Y,qX,qY,zN);

	ret = true;
	
	read_finish:
	mpq_clrs(qX,qY);
	
	return ret;
}
