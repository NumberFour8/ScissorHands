#include "generators.h"
#include "helpers.h"

FileGenerator::FileGenerator(mpz_t n,string filename)
 : BasicGenerator(n)
{
	fp.open(file);
}

FileGenerator::~FileGenerator()
{
	if (fp.is_open()) fp.close();
}

bool FileGenerator::next(ReducedPoint& P)
{	
	string ln;
	bool ret = true;
	
	mpq_t qX,qY;
	mpq_intz(qX,qY);
	
	if (!fp.is_open() || !getline(fp,ln)) 
	{
		cout << "ERROR: Cannot read from file." << endl;
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

	A = (ln == "-1" ? -1 : 1);

	// Precti racionalni X-ovou souradnici 
	fp >> ln;
	mpq_set_str(qX,ln.c_str(),10);
	
	// Precti racionalni Y-ovou souradnici 
	fp >> ln;
	mpq_set_str(qY,ln.c_str(),10);

	// Pokus se X-ovou a Y-ovou souradnici rekudovat modulo N
	try 
	{
		reduce_rational_point(P.X,P.Y,qX,qY,N));
		s++;
	}
	catch (mpz_t f)
	{
		cout << "ERROR: Cannot reduce on curve #" << s << endl;
		if (mpz_cmp_ui(f,0) != 0) // Byl nalezen faktor?
		{
			cout << "Factor found: " << mpz_to_string(f) << endl;
		}
		ret = false;
		currentA = 0;
		goto read_finish;
	}

	read_finish:
	mpq_clrs(qX,qY);
	
	return ret;
}
