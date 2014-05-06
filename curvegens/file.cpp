#include "generators.h"
#include "helpers.h"

FileGenerator::FileGenerator(mpz_t n,string filename)
 : BasicGenerator(n)
{
	fp.open(file);
	if (!fp.is_open())
	{
		cout << "ERROR: Cannot open file." << endl;
	} 
}

FileGenerator::~FileGenerator()
{
	if (fp.is_open()) fp.close();
}

int FileGenerator::getCurrentA()
{
	return currentA;
}

bool FileGenerator::next(ReducedPoint& P)
{
	string ln;
	bool ret = true;
	
	mpq_t qX,qY;
	mpq_intz(qX,qY);
	
	if (!getline(fp,ln)) 
	{
		ret = false;
		goto read_finish;
	}
	
	// Preskoc segment, ktery nezacina #
	if (ln.find("#") == string::npos) continue;
	
	// Je to prekroucena Edwardsova krivka s a = -1 ? 
	fp >> ln;
	if (ln != "-1" && ln != "1")
	{
		cout << "ERROR: Unsupported curve type." << endl;
		ret = false;
		goto read_finish;
	}

	currentA = (ln == "-1" ? -1 : 1);

	// Precti racionalni X-ovou souradnici 
	fp >> ln;
	mpq_set_str(qX,ln.c_str(),10);
	
	// Precti racionalni Y-ovou souradnici 
	fp >> ln;
	mpq_set_str(qY,ln.c_str(),10);

	// Pokus se X-ovou a Y-ovou souradnici rekudovat modulo N
	if (!reduce_rational_point(P.X,P.Y,qX,qY,N))
	{
		cout << "ERROR: Cannot reduce on curve #" << s << endl;
		if (mpz_cmp_ui(P.X,0) != 0) // Byl nalezen faktor?
		{
			cout << "Factor found: " << mpz_to_string(P.X) << endl;
		}
		else if (mpz_cmp_ui(P.Y,0) != 0)
		{
			cout << "Factor found: " << mpz_to_string(P.Y) << endl;
		}
		ret = false;
		currentA = 0;
		goto read_finish;
	}
	s++;

	read_finish:
	mpq_clrs(qX,qY);
	
	return ret;
}
