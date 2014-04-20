#include "helpers.h"

void to_mont_repr (mpz_t x, mpz_t n)
{
  mpz_mul_2exp (x, x, MAX_BITS);
  mpz_mod(x, x, n);
}

void from_mont_repr (mpz_t x, mpz_t n,mpz_t invB)
{
  mpz_mul(x, x, invB);
  mpz_mod(x, x, n);
}

void mpz_to_biguint (biguint_t a, mpz_t b)
{
  for (int i = 0;i < NB_DIGITS;i++)
  {
	#if GMP_NUMB_BITS == 32
		a[i] = mpz_getlimbn(b, i);  
	#else 
		if (i%2 == 0)
		  a[i] = (mpz_getlimbn(b, i/2) & 0x00000000ffffffff);
		else
		  a[i] = (mpz_getlimbn(b, i/2) >> 32);  
	#endif
  }
}

void biguint_to_mpz (mpz_t a, biguint_t b)
{
  int i;
  
  mpz_set_ui(a, 0);

  for (i = NB_DIGITS-1;i>=0;i--)
  {
    mpz_mul_2exp(a, a, 32);
	mpz_add_ui(a, a, b[i]);
  }
}

string mpz_to_string(mpz_t number)
{
	size_t sz = mpz_sizeinbase(number,10)+3;
	char* out = new char[sz];
	memset((void*)out,0,sz);
	
	mpz_get_str((char*)out,10,number);
	
	string ret(out);
	delete[] out;

	return ret;
}

void printBigInt(const char* tag,biguint_t B)
{
	printf("%s: {",tag);
	for (int i = 0;i < NB_DIGITS;++i)
	{
		printf("%#010x",B[i]);
		if (B[i+1] == 0) break;
		if (i != NB_DIGITS-1) printf(",");
	}
	printf("}\n");
}

void lcmToN(mpz_t res,const unsigned int n)
{
	mpz_set_ui(res,1);
	for (unsigned int i = 2;i <= n;++i)
	{
		mpz_lcm_ui(res,res,i);
	}
}

bool try_invert_mod(mpz_t invx,mpz_t x,mpz_t N)
{
	bool ret = false;
		 
	mpz_t b,r;
	mpz_intz(b,r);
         
	mpz_gcdext(r,invx,b,x,N);
    if (mpz_cmp_ui(r,1) == 0)
	{
		if (mpz_sgn(invx) == -1) 
		  mpz_add(invx,invx,N);
		ret = true;
	}
    else if (mpz_cmp(r,N) == 0)
	{
		// Pri chybe do r nastav 0
		mpz_set_si(invx,0);
	}
	else 
	{
		// Do r nastav faktor cisla N
	    mpz_set(invx,r); 
	}

	mpz_clrs(r,b);
	
	return ret;
}

bool reduce_mod(mpz_t r,mpq_t q,mpz_t n)
{
	mpz_t den,num;
	mpz_init_set(den,mpq_denref(q));
	mpz_init_set(num,mpq_numref(q));
	
	mpz_mod(den,den,n);
	mpz_mod(num,num,n);

	if (mpz_sgn(den) == -1) mpz_add(den,den,n);
	if (mpz_sgn(num) == -1) mpz_add(num,num,n);

	bool s = try_invert_mod(r,den,n);
	mpz_clear(den);

	if (!s) 
	{
	  // Pri selhani je v r zapsan faktor cisla N nebo 0, viz try_invert_mod().
	  mpz_clear(num);
	  return false; 
	}
	
	mpz_mul(r,r,num);
	mpz_mod(r,r,n);

	mpz_clear(num);
	return true;
}

void NAF::initialize(mpz_t N,unsigned char W)
{
	mpz_t number;
	mpz_init_set(number,N);
	w = W;
	
	if (bits != NULL) free((void*)bits);

	size_t sz = mpz_sizeinbase(number,2);
	bits = (char*)malloc(sz*8+1);
	l = 0;

	unsigned int pow2 = 1 << w,rem,pww = 1 << (w-1);

	for (int i = 0;mpz_sgn(number) > 0;++i)
	{
		if (mpz_odd_p(number))
		{
			rem = (unsigned int)mpz_fdiv_ui(number,pow2);
			if (rem&pww)
			{
				bits[i] = (unsigned char)(rem - pow2);
			}
			else {
				bits[i] = (unsigned char)rem;
			}
			if (bits[i] < 0)
			  mpz_add_ui(number,number,(int)(bits[i]*(-1)));
			else mpz_sub_ui(number,number,(int)(bits[i]));
			l = i;
		}
		else bits[i] = 0;
		mpz_div_2exp(number,number,1);
	}
	mpz_clear(number);

	l += 1;
	bits = (char*)realloc((void*)bits,(size_t)l*8);
}

void NAF::print() const
{
	if (bits == NULL)
	{
		printf("NAF not initialized");
		return;
	}

	printf("NAF%d: ",(int)w);
	for (int i = l-1;i >= 0;--i)
	{
		printf("%d",(int)bits[i]);
	}
	printf("\n");
}

NAF::~NAF()
{
	if (bits != NULL)
	{ 
	  free((void*)bits);
	  bits = NULL;
	}
}

// Precte pocátecni body krivek v afinnich souradnicich ze souboru
int readCurves(string file,mpz_t N,ExtendedPoint** pInit,bool& minus1)
{
	ifstream fp;
	if (file.length() < 2) 
	{
		file = "curves_twisted.txt";
		cout << "INFO: Defaulting to " << file  << endl;	
	}
	
	// Pokud se nepodari otevrit soubor, skonci
	fp.open(file);
	if (!fp.is_open())
	{
		cout << "ERROR: Cannot open file." << endl;
		return 0;
	} 
	
	// Inicializace racionalnich souradnic
	mpq_t qX,qY;
	mpq_intz(qX,qY);

	// Inicializace celociselnych souradnic
	mpz_t zX,zY;
	mpz_init_set_ui(zX,0);
	mpz_init_set_ui(zY,0);
	
	string ln;
	vector<ExtendedPoint> v;
	minus1 = true;
	cout << "Loading curves..." << endl;
	while (getline(fp,ln))
	{
		// Preskoc segment, ktery nezacina #
		if (ln.find("#") == string::npos) continue;
		
		// Je to prekroucena Edwardsova krivka s a = -1 ? 
		fp >> ln;
		if (v.size() > 0 && (ln == "-1") != minus1)
		{
			cout << "ERROR: Cannot mix curves with a = 1 and curves with a = -1." << endl; 
		
			fp.close();
			mpz_clrs(zX,zY);
			mpq_clrs(qX,qY);

			return 0;
		}
		minus1 = (ln == "-1");
		
		// Precti racionalni X-ovou souradnici a zkrat
		fp >> ln;
		mpq_set_str(qX,ln.c_str(),10);
		mpq_canonicalize(qX);
		
		// Precti racionalni Y-ovou souradnici a zkrat
		fp >> ln;
		mpq_set_str(qY,ln.c_str(),10);
		mpq_canonicalize(qY);

		// Pokus se X-ovou a Y-ovou souradnici rekudovat modulo N
		if (!reduce_mod(zX,qX,N) || !reduce_mod(zY,qY,N))
		{
			cout << "ERROR: Cannot reduce on curve #" << v.size() << endl;
			if (mpz_cmp_ui(zX,0) != 0) // Byl nalezen faktor?
			{
				cout << "Factor found: " << mpz_to_string(zX) << endl;
			}
			else if (mpz_cmp_ui(zY,0) != 0)
			{
				cout << "Factor found: " << mpz_to_string(zY) << endl;
			}
			
			fp.close();
			mpz_clrs(zX,zY);
			mpq_clrs(qX,qY);

			return 0;
		}

		// Vytvor bod v Extended souradnicích z redukovanych afinnich bodu modulo N
		v.push_back(ExtendedPoint(zX,zY,N)); 
	}

	// Prekopiruj body z vektoru do pameti
	*pInit = new ExtendedPoint[v.size()];
	std::copy(v.begin(),v.end(),*pInit);

	// Vycisti pamet
	fp.close();
	mpz_clrs(zX,zY);
	mpq_clrs(qX,qY);

	return (int)v.size();
}
