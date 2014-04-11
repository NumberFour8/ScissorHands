#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

#include "kernel.h"

using namespace std;

int main()
{
	string cf,X,Y,N;
	mpz_t zcf,zX,zY,zN;

	ifstream fp;
	fp.open("C:\\Users\\Luc\\Desktop\\edwards.txt");

	fp >> N;
	cout << "N = " << N << endl;
	fp >> X;
	cout << "X = " << X << endl;
	fp >> Y;
	cout << "Y = " << Y << endl;
	fp >> cf;
	cout << "B = " << cf << endl;
	cout << "a = -1" << endl << endl;
	fp.close();
	
	cout << "Is this okay ?" << endl;
	char c;
	cin >> c;
	if (c == 'n') return 0;
	
	////////////////////////////////////////////////////////////
	
	mpz_init_set_str(zN,N.c_str(),10);
	mpz_init_set_str(zX,X.c_str(),10);
	mpz_init_set_str(zY,Y.c_str(),10);
	
	// Koefcient v NAF rozvoji, do kterého se chceme dopočítat
	mpz_init(zcf);
	//mpz_set_ui(zcf,(unsigned int)std::stoul(cf));
	lcmToN(zcf,(unsigned int)std::stoul(cf));

	cout << "lcm(1,2,..." << cf << ") = ";
	printmpz("%s\n",zcf);

	NAF coeffNaf(2,zcf);
	//coeffNaf.print();
	mpz_clear(zcf);
		
	ExtendedPoint pts;
	pts.fromAffine(zX,zY,zN);

	ExtendedPoint neutral;
	neutral.infinity(zN);

	// Spočítáme W = 2^32, 3*N, -N^(-1) mod W a W^(-1) mod N
	mpz_t z3N,zInvW,zInvN;
	mpz_init(z3N);
	mpz_mul_ui(z3N,zN,3);
	mpz_init(zInvW);
	mpz_ui_pow_ui (zInvW, 2, SIZE_DIGIT); 
    mpz_init(zInvN);	
	mpz_invert (zInvN, zN, zInvW);
    mpz_sub (zInvN, zInvW, zInvN);
	mpz_invert (zInvW, zInvW, zN);

	// Pomocná struktura
	Aux ax;
	memset(&ax,0,sizeof(Aux));

	mpz_to_biguint(ax.N,zN);
	mpz_to_biguint(ax.N3,z3N);
	ax.invN = (digit_t)mpz_get_ui(zInvN);

	cudaError_t cudaStatus = compute(ax,&neutral,&pts,coeffNaf);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

	cout << endl;

	
	ExtResult res;
	if (pts.toAffine(zX,zY,zN,zInvW,&res)){
		printmpz("(%s,",zX);
		printmpz("%s)\n",zY);
	}
	else {
		printmpz("Factor: %s",res.factor);
	}
    
	mpz_clear(zX);
	mpz_clear(zY);
	
	mpz_clear(zN);
	mpz_clear(z3N);
	mpz_clear(zInvW);
	mpz_clear(zInvN);

	cin >> c;
    return 0;
}
