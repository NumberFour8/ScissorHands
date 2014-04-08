#include <iostream>
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
	fp >> k;
	cout << "k = " << k << endl;
	cout << "a = -1" << endl;
	fp.close();
	
	cout << "Is this okay ?" << endl;
	char c;
	cin >> c;
	if (c == 'n') return 0;
	
	////////////////////////////////////////////////////////////
	
	mpz_init_set_str(zN,N.c_str(),10);
	mpz_init_set_str(zcf,cf.c_str(),10);
	mpz_init_set_str(zX,X.c_str(),10);
	mpz_init_set_str(zY,Y.c_str(),10);
	
	// Koefcient v NAF rozvoji, do kterého se chceme dopočítat
	NAF coeffNaf(2,zcf);
	mpz_clear(zcf);
		
	h_ExtendedPoint pts;
	pts.fromAffine(zX,zY,zN);

	// Spočítáme W = 2^32, 3*N, -N^(-1) mod W a W^(-1) mod N
	mpz_t z3N,zInvW,zInvN;
	mpz_inits(z3N,zInvW,zInvN);
	mpz_mul_ui(z3N,zN,3);
	mpz_ui_pow_ui (zInvW, 2, SIZE_DIGIT); 
    mpz_invert (zInvN, zN, zInvW);
    mpz_sub (zInvN, zInvW, zInvN);
	mpz_invert (zInvW, zInvW, zN);

	// Pomocná struktura
	h_Aux ax;
	mpz_to_biguint(ax.N,zN);
	mpz_to_biguint(ax.N3,z3N);
	ax.invN = (digit_t)mpz_get_ui(zInvN);

	cudaError_t cudaStatus = computeExtended(ax,&pts,coeffNaf);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
    
	mpz_clears(zX,zY,zN);
	mpz_clears(z3N,zInvW,zInvN);

	cin >> c;
    return 0;
}
